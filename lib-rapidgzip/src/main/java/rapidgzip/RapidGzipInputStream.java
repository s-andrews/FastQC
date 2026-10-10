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
 *
 * Pipeline: {@link BlockFinder} guesses a block start in each chunk,
 * {@link DeflateDecoder} decodes from it into a {@link Chunk} (16 bit markers
 * while the unknown window may still be referenced), and this class orders the
 * chunks, resolves their markers and checks the gzip trailers.
 *
 * Why the output is identical to a sequential decoder: decoding from a real
 * block boundary with the right window is deterministic, and a chunk is only
 * accepted when its start bit equals the end bit reached by decoding the
 * previous output, i.e. when its start is proven to be a real boundary. A false
 * block start either fails to decode (DataFormatError/EOF, so the finder tries
 * the next candidate) or produces a chunk whose start doesn't match, which is
 * dropped and that range decoded sequentially. The CRC32 and ISIZE of each
 * member are checked on top of that.
 *
 * Threading: one caller thread reads, orders chunks, decodes sequentially when
 * needed and checks CRCs; a fixed pool of {@code threads} daemon workers finds
 * blocks, decodes chunks and resolves markers. Only {@link #ratio},
 * {@link #expectedBit} and the {@link Buffers} pool are shared with other
 * threads. Workers stop at the end of the data, or on {@link #close()}.
 *
 * Memory: in FASTQ ~55% of a speculatively decoded chunk refers back to its
 * unknown window, so chunks usually stay in marker mode (2 bytes per output
 * byte) until resolved: ~20 bytes per compressed byte. Chunks being decoded
 * plus chunks waiting to be read are capped at a quarter of the max heap.
 *
 * Benefit: a 510MB AVITI file decompresses in 1.89s on 4 threads and 1.42s on
 * 8, vs 3.58s with GZIPInputStream, at ~2.3x zlib's total CPU.
 *
 * Not thread safe, like any InputStream; needs random access, so not for pipes.
 */
public final class RapidGzipInputStream extends InputStream {

	// Compressed bytes per chunk: small enough to bound memory in marker mode (see
	// the constructor), large enough that block finding and the sequential
	// fallbacks are a small part of the work.
	public static final int DEFAULT_CHUNK_SIZE = 1 << 20;

	// Bytes read from the file at a time by each BitReader
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
	// Heap for the output of the chunks in flight
	private final long budget;

	// Chunk output buffers, reused across the decoders, the resolving tasks and the reader
	private final Buffers buffers;
	// Speculative decodes of chunks [nextChunk, nextToSubmit), in chunk order: the
	// head is always the result for chunk nextChunk
	private final ArrayDeque<Future<Chunk>> inFlight = new ArrayDeque<Future<Chunk>>();
	private int nextToSubmit = 0;
	// Index of the next chunk to put in order
	private int nextChunk = 0;
	private final long firstBlockBit;
	// Expected decompressed bytes per compressed byte, to size output buffers
	private final AtomicInteger ratio = new AtomicInteger(4);

	// Bit offset of the next block boundary to decode, and the output preceding it.
	// Volatile only so that compressedPosition() can be called from another thread.
	private volatile long expectedBit;
	private byte [] window = new byte[0];
	private boolean finished = false;
	// The first failure reading the data: every later read fails with it too, rather
	// than returning data past a corrupt member
	private IOException failure = null;

	// CRC32 and size of the output of the current member so far
	private final CRC32 crc = new CRC32();
	private long memberSize = 0;

	// Output ready to be read, as ranges of chunk arrays (no copy)
	private final ArrayDeque<byte []> segments = new ArrayDeque<byte []>();
	private final ArrayDeque<int []> segmentRanges = new ArrayDeque<int []>();
	private byte [] current = null;
	private int currentPos;
	private int currentEnd;
	private boolean closed = false;

	/**
	 * Opens a gzip file, checks its first header and starts decoding the first
	 * chunks on {@code threads} worker threads.
	 *
	 * @throws ZipException if the file doesn't start with a gzip header
	 * @throws EOFException if the file is empty
	 */
	public RapidGzipInputStream (Path file, int threads) throws IOException {
		this(open(file, threads, DEFAULT_CHUNK_SIZE), true, threads, DEFAULT_CHUNK_SIZE);
	}

	public RapidGzipInputStream (Path file, int threads, int chunkSize) throws IOException {
		this(open(file, threads, chunkSize), true, threads, chunkSize);
	}

	// Checks the arguments before opening the file, so that bad ones don't leak it
	private static FileChannel open (Path file, int threads, int chunkSize) throws IOException {
		if (threads < 1) throw new IllegalArgumentException("threads must be >= 1");
		if (chunkSize < 64 << 10) throw new IllegalArgumentException("chunkSize must be >= 64KB");
		return FileChannel.open(file, StandardOpenOption.READ);
	}

	private RapidGzipInputStream (FileChannel channel, boolean closeChannel, int threads, int chunkSize) throws IOException {
		this.channel = channel;
		this.closeChannel = closeChannel;
		try {
			this.fileSize = channel.size();
		}
		catch (IOException e) {
			if (closeChannel) channel.close();
			throw e;
		}
		this.chunkSize = chunkSize;
		this.chunkCount = (int) Math.max(1, (fileSize + chunkSize - 1) / chunkSize);
		this.source = this::readInto;
		// Each chunk in flight holds its output, often all as 16 bit values: for FASTQ ~6x the chunk
		// size as 16 bit values, then as bytes. Use up to a quarter of the heap for them.
		// (A 2MB chunk, 2 per thread, ran a 512MB heap out of memory.) The 20x estimate assumes a
		// FASTQ-like compression ratio: much more compressible data uses proportionally more.
		long perChunk = (long) chunkSize * 20;
		this.budget = Runtime.getRuntime().maxMemory() / 4;
		this.maxInFlight = (int) Math.max(2, Math.min(threads + 2, budget / perChunk));
		this.buffers = new Buffers(maxInFlight);

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

		// Daemon threads, so an unclosed stream doesn't keep the JVM alive
		pool = Executors.newFixedThreadPool(threads, r -> {
			Thread t = new Thread(r, "rapidgzip");
			t.setDaemon(true);
			return t;
		});
		submitChunks();
	}

	/**
	 * Positional read into dest, safe to call from several threads at once as it
	 * doesn't use the channel position. Returns the number of bytes read, fewer
	 * than length only at the end of the file.
	 */
	private int readInto (long offset, byte [] dest, int destOffset, int length) throws IOException {
		int n = (int) Math.max(0, Math.min(length, fileSize - offset));
		ByteBuffer buffer = ByteBuffer.wrap(dest, destOffset, n);
		while (buffer.hasRemaining()) {
			if (channel.read(buffer, offset + (buffer.position() - destOffset)) < 0) break;
		}
		return buffer.position() - destOffset;
	}

	// Compressed data of the chunk a worker thread is looking for a block start in
	private final ThreadLocal<byte []> chunkData = new ThreadLocal<byte []>();

	/** Bit where chunk k's range starts, which is where chunk k-1 stops decoding (at the next boundary). */
	private long chunkStartBit (int k) {
		return k >= chunkCount ? Long.MAX_VALUE : (long) k * chunkSize * 8;
	}

	// Output buffer size for a chunk from the ratio seen so far, plus 1/8, so that it's
	// rarely grown (a copy of the largest array in use)
	private int sizeHint () {
		return (int) Math.min(Integer.MAX_VALUE - 64, (long) chunkSize * ratio.get() * 9 / 8);
	}

	/**
	 * Chunks allowed in flight: maxInFlight, fewer if the data compresses better than
	 * FASTQ. In marker mode a chunk holds ~3 bytes per output byte (16 bit values, then
	 * resolved bytes), so from the ratio seen so far each takes chunkSize * 3 * ratio,
	 * against the chunkSize * 20 maxInFlight assumes. Unchanged for FASTQ (ratio <= 6).
	 */
	private int inFlightLimit () {
		long perChunk = (long) chunkSize * Math.max(20, 3L * ratio.get());
		return (int) Math.max(2, Math.min(maxInFlight, budget / perChunk));
	}

	// Chunks being decoded and chunks decoded but not yet read both hold their output,
	// so both count against the in-flight limit
	private void submitChunks () {
		while (nextToSubmit < chunkCount && inFlight.size() + ordered.size() < inFlightLimit()) {
			final int k = nextToSubmit++;
			inFlight.add(pool.submit(() -> decodeChunk(k)));
		}
	}

	/**
	 * Speculatively decodes chunk k, or returns null if no block boundary was found.
	 * Runs on a worker thread. Candidates from the finder are tried in order until
	 * one decodes without error up to the next chunk; a false positive almost
	 * always fails within a few symbols. A candidate that decodes but isn't the
	 * real boundary is caught later by the start bit check in orderNext().
	 */
	private Chunk decodeChunk (int k) throws IOException {
		long stopBit = chunkStartBit(k+1);
		if (k == 0) {
			// The first block start and its (empty) window are known: no markers
			Chunk chunk = new DeflateDecoder(new BitReader(source, SEGMENT_SIZE), buffers).decode(firstBlockBit, stopBit, new byte[0], sizeHint());
			updateRatio(chunk);
			return chunk;
		}

		long start = (long) k * chunkSize;
		int length = (int) Math.min(chunkSize, fileSize - start);
		// Candidates are searched in this chunk only; the slack lets headers starting near its end be checked
		byte [] data = chunkData.get();
		if (data == null || data.length < chunkSize + HEADER_SLACK) {
			data = new byte[chunkSize + HEADER_SLACK];
			chunkData.set(data);
		}
		int read = readInto(start, data, 0, length + HEADER_SLACK);
		BlockFinder finder = new BlockFinder();
		long from = 0;
		while (true) {
			long found = finder.find(data, read, from, (long) length * 8);
			if (found < 0) return null;
			try {
				Chunk chunk = new DeflateDecoder(new BitReader(source, SEGMENT_SIZE), buffers).decode(start * 8 + found, stopBit, null, sizeHint());
				updateRatio(chunk);
				return chunk;
			}
			catch (DataFormatError | EOFException e) {
				// Not a real block boundary
				from = found + 1;
			}
		}
	}

	// Last ratio seen by any worker (a benign race), capped at 64 for the buffer sizes
	private void updateRatio (Chunk chunk) {
		long compressed = Math.max(1, (chunk.endBit - chunk.startBit) / 8);
		ratio.set((int) Math.max(1, Math.min(64, (chunk.length() + compressed - 1) / compressed)));
	}

	/**
	 * Waits for the speculative result of chunk nextChunk, or null if it failed.
	 * A worker exception other than an Error isn't reported here: the range is
	 * decoded sequentially, which reports a real problem with the data.
	 */
	private Chunk nextResult () throws IOException {
		if (inFlight.isEmpty()) submitChunks();
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

	// On the caller thread, from the known boundary with the known window: at least one block, then up
	// to the first boundary at or after stopBit
	private Chunk decodeSequentially (long stopBit) throws IOException {
		sequentialDecodes++;
		return new DeflateDecoder(new BitReader(source, SEGMENT_SIZE), buffers).decode(expectedBit, stopBit, window, sizeHint());
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

	// Chunks put in order, waiting for their markers to be resolved and to be read
	private final ArrayDeque<Ordered> ordered = new ArrayDeque<Ordered>();

	// For tests: chunks used as decoded speculatively, and sequential decodes when speculation failed
	int speculativeChunks = 0;
	int sequentialDecodes = 0;

	/**
	 * Takes the next chunk in order: a speculative one if it starts where the
	 * previous one ended, else decodes sequentially. Its window gives the window
	 * of the next one, and its markers are resolved on the pool.
	 */
	private void orderNext () throws IOException {
		if (nextChunk >= chunkCount) {
			// All chunks used or dropped but the data hasn't ended: decode the rest
			addOrdered(decodeSequentially(Long.MAX_VALUE));
			return;
		}

		Chunk chunk = nextResult();
		// The correctness check: decoding the known stream reached exactly this chunk's
		// start, so it's a real block boundary and its output is what a sequential decoder gives
		if (chunk != null && chunk.startBit == expectedBit) {
			nextChunk++;
			if (nextChunk > 1) speculativeChunks++;
			addOrdered(chunk);
			return;
		}

		// Speculation failed: decode sequentially up to its start, or else over the whole chunk.
		// A start after expectedBit is usually a real boundary the previous chunk's decode
		// stopped short of (the first boundary was a stored or fixed block, which the finder
		// doesn't look for). A start before it was a false positive. Each pass decodes at
		// least one block, so this always makes progress.
		boolean usable = chunk != null && chunk.startBit > expectedBit;
		Chunk sequential = decodeSequentially(usable ? chunk.startBit : chunkStartBit(nextChunk+1));
		addOrdered(sequential);
		if (usable && sequential.endBit == chunk.startBit && !finished) {
			// It matches now: put it back as the head result for chunk nextChunk
			inFlight.addFirst(java.util.concurrent.CompletableFuture.completedFuture(chunk));
		}
		else {
			nextChunk++;
		}
	}

	/**
	 * Appends a chunk to the output order. Only the last 32KB are resolved here, on
	 * the caller thread, as the next chunk needs them as its window; resolving the
	 * whole marker part goes to the pool so the caller can move on to the next chunk.
	 */
	private void addOrdered (Chunk chunk) throws IOException {
		final byte [] before = window;
		window = chunk.windowAfter(before);
		expectedBit = chunk.endBit;
		if (chunk.endOfData) finished = true;
		Future<byte []> resolved = chunk.markersTo > chunk.markersFrom
				? pool.submit(() -> chunk.resolveAndRelease(before, buffers))
				: java.util.concurrent.CompletableFuture.completedFuture(new byte[0]);
		ordered.add(new Ordered(chunk, resolved));
	}

	/** Queues the output of the next chunk(s), returns false at the end. */
	private boolean fill () throws IOException {
		while (segments.isEmpty()) {
			// Keep the pipeline busy, waiting for a speculative chunk only if there's nothing else to do
			while (!finished && ordered.size() < inFlightLimit()
					&& (ordered.isEmpty() || nextChunk >= chunkCount || (inFlight.peek() != null && inFlight.peek().isDone()))) {
				orderNext();
			}
			Ordered next = ordered.poll();
			if (next == null) return false;
			submitChunks();
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

	// A chunk's output is its resolved marker part followed by its byte part.
	// The resolved buffer comes from the pool, so may be longer than the markers
	private void emit (Chunk chunk, byte [] resolved) throws IOException {
		addSegment(resolved, 0, chunk.markersLength());
		addSegment(chunk.bytes, chunk.bytesFrom, chunk.bytesTo);
		checkMembers(chunk, resolved);
	}

	private void addSegment (byte [] data, int from, int to) {
		if (to > from) {
			segments.add(data);
			segmentRanges.add(new int [] { from, to });
		}
	}

	/**
	 * Updates the CRC32 and size of the current member and checks them at each member end.
	 * Member ends are output offsets within the chunk, which may fall in either part.
	 * Runs on the caller thread, in output order, as CRC32 can't be combined here.
	 */
	private void checkMembers (Chunk chunk, byte [] resolved) throws IOException {
		int m = 0;
		long offset = 0;
		byte [][] parts = { resolved, chunk.bytes };
		int [] froms = { 0, chunk.bytesFrom };
		int [] tos = { chunk.markersLength(), chunk.bytesTo };
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
		if (failure != null) throw failure;
		// Each buffer is in one segment, so it can be reused once read: its CRC was
		// checked in emit, and the window after it taken in addOrdered
		if (current != null) {
			buffers.release(current);
			current = null;
		}
		boolean more;
		try {
			more = fill();
		}
		catch (IOException e) {
			failure = e;
			pool.shutdownNow();
			throw e;
		}
		if (!more) {
			// The end of the data: stop the workers now rather than at close(), so a
			// stream that is read to the end but never closed doesn't keep them, nor
			// keep decoding chunks past trailing non-gzip data
			pool.shutdownNow();
			return false;
		}
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

	/**
	 * How far decompression has got through the compressed file, in bytes, for progress reporting.
	 * Chunks are ordered ahead of the reader, so this can lead what has been read by a few chunks.
	 * Can be called from any thread.
	 */
	public long compressedPosition () {
		return Math.min(fileSize, expectedBit >>> 3);
	}

	/** Size of the compressed file in bytes. */
	public long compressedSize () {
		return fileSize;
	}

	@Override
	public int available () {
		return current == null ? 0 : currentEnd - currentPos;
	}

	/**
	 * Stops the workers and closes the file. Doesn't wait for them: a worker still
	 * decoding fails on the closed (or, being interrupted, closed-by-interrupt)
	 * channel and its result is dropped.
	 */
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
