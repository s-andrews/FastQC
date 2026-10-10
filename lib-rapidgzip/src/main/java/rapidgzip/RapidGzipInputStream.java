package rapidgzip;

import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.ByteBuffer;
import java.nio.channels.FileChannel;
import java.nio.file.Path;
import java.nio.file.StandardOpenOption;
import java.util.ArrayDeque;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.atomic.AtomicInteger;
import java.util.zip.CRC32;
import java.util.zip.ZipException;

/**
 * Decompresses a gzip file using several threads, in the manner of
 * <a href="https://github.com/mxmlnkn/rapidgzip">rapidgzip</a>.
 *
 * The compressed file is split into chunks. For each chunk a worker thread
 * looks for the first deflate block boundary in it and decodes from there up to
 * the first block boundary in the next chunk, without knowing the 32KB of
 * output preceding it: references to that window are kept as markers. The
 * reader then takes the chunks in order: a chunk is only used if it starts
 * exactly where the previous one ended, which proves its start was a real block
 * boundary, and its markers are replaced with the now known window. Anything
 * which can't be matched that way is decoded sequentially. The CRC32 and size
 * of every gzip member are checked, so the output is exactly what a sequential
 * decoder produces.
 *
 * Multiple members are decoded as one stream, and anything after the last
 * member which isn't a gzip header is ignored, as GZIPInputStream does.
 */
public final class RapidGzipInputStream extends InputStream {

	public static final int DEFAULT_CHUNK_SIZE = 2 << 20;

	private static final int SEGMENT_SIZE = 256 << 10;
	// Max size of a block header, so the finder can check headers near the end of a chunk
	private static final int HEADER_SLACK = 1024;

	private final FileChannel channel;
	private final boolean closeChannel;
	private final long fileSize;
	private final int chunkSize;
	private final int chunkCount;
	private final BitReader.Source source;
	private final ExecutorService pool;
	private final int maxInFlight;

	private final ArrayDeque<Future<Chunk>> inFlight = new ArrayDeque<Future<Chunk>>();
	private int nextToSubmit = 0;
	private int nextChunk = 0;
	private final long firstBlockBit;
	// Expected decompressed bytes per compressed byte, to size output buffers
	private final AtomicInteger ratio = new AtomicInteger(4);

	// Bit offset of the next block boundary to decode, and the output preceding it
	private long expectedBit;
	private byte [] window = new byte[0];
	private boolean finished = false;

	private final CRC32 crc = new CRC32();
	private long memberSize = 0;

	private final ArrayDeque<byte []> segments = new ArrayDeque<byte []>();
	private final ArrayDeque<int []> segmentRanges = new ArrayDeque<int []>();
	private byte [] current = null;
	private int currentPos;
	private int currentEnd;
	private boolean closed = false;

	public RapidGzipInputStream (Path file, int threads) throws IOException {
		this(FileChannel.open(file, StandardOpenOption.READ), true, threads, DEFAULT_CHUNK_SIZE);
	}

	public RapidGzipInputStream (Path file, int threads, int chunkSize) throws IOException {
		this(FileChannel.open(file, StandardOpenOption.READ), true, threads, chunkSize);
	}

	private RapidGzipInputStream (FileChannel channel, boolean closeChannel, int threads, int chunkSize) throws IOException {
		if (threads < 1) throw new IllegalArgumentException("threads must be >= 1");
		if (chunkSize < 64 << 10) throw new IllegalArgumentException("chunkSize must be >= 64KB");
		this.channel = channel;
		this.closeChannel = closeChannel;
		this.fileSize = channel.size();
		this.chunkSize = chunkSize;
		this.chunkCount = (int) Math.max(1, (fileSize + chunkSize - 1) / chunkSize);
		this.source = this::readAt;
		this.maxInFlight = threads * 2;

		try {
			if (fileSize == 0) throw new EOFException();
			BitReader header = new BitReader(source, HEADER_SLACK);
			header.seek(0);
			if (!GzipHeader.skip(header)) throw new ZipException("Not in GZIP format");
			firstBlockBit = header.position();
		}
		catch (IOException e) {
			if (closeChannel) channel.close();
			throw e;
		}
		expectedBit = firstBlockBit;

		pool = Executors.newFixedThreadPool(threads, r -> {
			Thread t = new Thread(r, "rapidgzip");
			t.setDaemon(true);
			return t;
		});
		submitChunks();
	}

	private byte [] readAt (long offset, int length) throws IOException {
		int n = (int) Math.max(0, Math.min(length, fileSize - offset));
		ByteBuffer buffer = ByteBuffer.allocate(n);
		while (buffer.hasRemaining()) {
			if (channel.read(buffer, offset + buffer.position()) < 0) break;
		}
		return buffer.position() == n ? buffer.array() : java.util.Arrays.copyOf(buffer.array(), buffer.position());
	}

	private long chunkStartBit (int k) {
		return k >= chunkCount ? Long.MAX_VALUE : (long) k * chunkSize * 8;
	}

	private int sizeHint () {
		return (int) Math.min(Integer.MAX_VALUE - 64, (long) chunkSize * ratio.get() * 9 / 8);
	}

	private void submitChunks () {
		while (nextToSubmit < chunkCount && inFlight.size() < maxInFlight) {
			final int k = nextToSubmit++;
			inFlight.add(pool.submit(() -> decodeChunk(k)));
		}
	}

	/** Speculatively decodes chunk k, or returns null if no block boundary was found. */
	private Chunk decodeChunk (int k) throws IOException {
		long stopBit = chunkStartBit(k+1);
		if (k == 0) {
			Chunk chunk = new DeflateDecoder(new BitReader(source, SEGMENT_SIZE)).decode(firstBlockBit, stopBit, new byte[0], sizeHint());
			updateRatio(chunk);
			return chunk;
		}

		long start = (long) k * chunkSize;
		int length = (int) Math.min(chunkSize, fileSize - start);
		byte [] data = readAt(start, length + HEADER_SLACK);
		BlockFinder finder = new BlockFinder();
		long from = 0;
		while (true) {
			long found = finder.find(data, data.length, from, (long) length * 8);
			if (found < 0) return null;
			try {
				Chunk chunk = new DeflateDecoder(new BitReader(source, SEGMENT_SIZE)).decode(start * 8 + found, stopBit, null, sizeHint());
				updateRatio(chunk);
				return chunk;
			}
			catch (DataFormatError | EOFException e) {
				// Not a real block boundary
				from = found + 1;
			}
		}
	}

	private void updateRatio (Chunk chunk) {
		long compressed = Math.max(1, (chunk.endBit - chunk.startBit) / 8);
		ratio.set((int) Math.max(1, Math.min(64, (chunk.length() + compressed - 1) / compressed)));
	}

	private Chunk nextResult () throws IOException {
		Future<Chunk> future = inFlight.poll();
		submitChunks();
		try {
			return future.get();
		}
		catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new java.io.InterruptedIOException();
		}
		catch (ExecutionException e) {
			if (e.getCause() instanceof Error) throw (Error) e.getCause();
			// Decoded sequentially instead, which reports any real error
			return null;
		}
	}

	private Chunk decodeSequentially (long stopBit) throws IOException {
		sequentialDecodes++;
		return new DeflateDecoder(new BitReader(source, SEGMENT_SIZE)).decode(expectedBit, stopBit, window, sizeHint());
	}

	/** A chunk in order, with its markers being resolved. */
	private static final class Ordered {
		final Chunk chunk;
		final Future<byte []> resolved;

		Ordered (Chunk chunk, Future<byte []> resolved) {
			this.chunk = chunk;
			this.resolved = resolved;
		}
	}

	private final ArrayDeque<Ordered> ordered = new ArrayDeque<Ordered>();

	// Chunks used as decoded speculatively, and sequential decodes when speculation failed
	int speculativeChunks = 0;
	int sequentialDecodes = 0;

	/**
	 * Takes the next chunk in order: a speculative one if it starts where the
	 * previous one ended, else decodes sequentially. Its window gives the window
	 * of the next one, and its markers are resolved on the pool.
	 */
	private void orderNext () throws IOException {
		if (nextChunk >= chunkCount) {
			addOrdered(decodeSequentially(Long.MAX_VALUE));
			return;
		}

		Chunk chunk = nextResult();
		if (chunk != null && chunk.startBit == expectedBit) {
			nextChunk++;
			if (nextChunk > 1) speculativeChunks++;
			addOrdered(chunk);
			return;
		}

		// Speculation failed: decode sequentially up to its start, or else over the whole chunk
		boolean usable = chunk != null && chunk.startBit > expectedBit;
		Chunk sequential = decodeSequentially(usable ? chunk.startBit : chunkStartBit(nextChunk+1));
		addOrdered(sequential);
		if (usable && sequential.endBit == chunk.startBit && !finished) {
			// It matches now
			inFlight.addFirst(java.util.concurrent.CompletableFuture.completedFuture(chunk));
		}
		else {
			nextChunk++;
		}
	}

	private void addOrdered (Chunk chunk) throws IOException {
		final byte [] before = window;
		window = chunk.windowAfter(before);
		expectedBit = chunk.endBit;
		if (chunk.endOfData) finished = true;
		Future<byte []> resolved = chunk.markersTo > chunk.markersFrom
				? pool.submit(() -> chunk.resolve(before))
				: java.util.concurrent.CompletableFuture.completedFuture(new byte[0]);
		ordered.add(new Ordered(chunk, resolved));
	}

	/** Queues the output of the next chunk(s), returns false at the end. */
	private boolean fill () throws IOException {
		while (segments.isEmpty()) {
			// Keep the pipeline busy, waiting for a speculative chunk only if there's nothing else to do
			while (!finished && ordered.size() < maxInFlight
					&& (ordered.isEmpty() || nextChunk >= chunkCount || (inFlight.peek() != null && inFlight.peek().isDone()))) {
				orderNext();
			}
			Ordered next = ordered.poll();
			if (next == null) return false;
			emit(next.chunk, get(next.resolved));
		}
		return true;
	}

	private byte [] get (Future<byte []> future) throws IOException {
		try {
			return future.get();
		}
		catch (InterruptedException e) {
			Thread.currentThread().interrupt();
			throw new java.io.InterruptedIOException();
		}
		catch (ExecutionException e) {
			if (e.getCause() instanceof IOException) throw (IOException) e.getCause();
			if (e.getCause() instanceof Error) throw (Error) e.getCause();
			throw new IOException(e.getCause());
		}
	}

	private void emit (Chunk chunk, byte [] resolved) throws IOException {
		addSegment(resolved, 0, resolved.length);
		addSegment(chunk.bytes, chunk.bytesFrom, chunk.bytesTo);
		checkMembers(chunk, resolved);
	}

	private void addSegment (byte [] data, int from, int to) {
		if (to > from) {
			segments.add(data);
			segmentRanges.add(new int [] { from, to });
		}
	}

	/** Updates the CRC32 and size of the current member and checks them at each member end. */
	private void checkMembers (Chunk chunk, byte [] resolved) throws IOException {
		int m = 0;
		long offset = 0;
		byte [][] parts = { resolved, chunk.bytes };
		int [] froms = { 0, chunk.bytesFrom };
		int [] tos = { resolved.length, chunk.bytesTo };
		for (int p = 0; p < 2; p++) {
			int pos = froms[p];
			int length = tos[p] - froms[p];
			while (m < chunk.memberEnds.length && chunk.memberEnds[m] <= offset + length) {
				int cut = froms[p] + (int) (chunk.memberEnds[m] - offset);
				crc.update(parts[p], pos, cut - pos);
				memberSize += cut - pos;
				endMember(chunk.memberCrcs[m], chunk.memberSizes[m]);
				pos = cut;
				m++;
			}
			crc.update(parts[p], pos, tos[p] - pos);
			memberSize += tos[p] - pos;
			offset += length;
		}
	}

	private void endMember (int expectedCrc, int expectedSize) throws ZipException {
		if ((int) crc.getValue() != expectedCrc || (int) memberSize != expectedSize) {
			throw new ZipException("Corrupt GZIP trailer");
		}
		crc.reset();
		memberSize = 0;
	}

	private boolean nextSegment () throws IOException {
		if (closed) throw new IOException("Stream closed");
		if (!fill()) return false;
		current = segments.poll();
		int [] range = segmentRanges.poll();
		currentPos = range[0];
		currentEnd = range[1];
		return true;
	}

	@Override
	public int read () throws IOException {
		if (current == null || currentPos == currentEnd) {
			if (!nextSegment()) return -1;
		}
		return current[currentPos++] & 0xFF;
	}

	@Override
	public int read (byte [] b, int off, int len) throws IOException {
		java.util.Objects.checkFromIndexSize(off, len, b.length);
		if (len == 0) return 0;
		if (current == null || currentPos == currentEnd) {
			if (!nextSegment()) return -1;
		}
		int n = Math.min(len, currentEnd - currentPos);
		System.arraycopy(current, currentPos, b, off, n);
		currentPos += n;
		return n;
	}

	@Override
	public int available () {
		return current == null ? 0 : currentEnd - currentPos;
	}

	@Override
	public void close () throws IOException {
		if (closed) return;
		closed = true;
		pool.shutdownNow();
		segments.clear();
		ordered.clear();
		segmentRanges.clear();
		current = null;
		if (closeChannel) channel.close();
	}

}
