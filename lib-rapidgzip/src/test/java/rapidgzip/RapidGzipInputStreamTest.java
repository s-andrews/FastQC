package rapidgzip;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertThrows;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.ByteArrayInputStream;
import java.io.ByteArrayOutputStream;
import java.io.EOFException;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Random;
import java.util.zip.CRC32;
import java.util.zip.Deflater;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.api.io.TempDir;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.CsvSource;
import org.junit.jupiter.params.provider.ValueSource;

public class RapidGzipInputStreamTest {

	private static final int MIN_CHUNK = 64 << 10;

	@TempDir
	Path dir;

	// --- Data ---

	/** Synthetic FASTQ: repetitive structure with random bases and qualities, like real reads. */
	static byte [] fastq (int reads, long seed) {
		Random random = new Random(seed);
		StringBuilder sb = new StringBuilder();
		String bases = "ACGT";
		for (int r = 0; r < reads; r++) {
			sb.append("@SRR000000.").append(r).append(" HWI-ST:8:1101:").append(random.nextInt(20000)).append(':').append(random.nextInt(200000)).append(" length=150\n");
			for (int i = 0; i < 150; i++) sb.append(random.nextInt(50) == 0 ? 'N' : bases.charAt(random.nextInt(4)));
			sb.append("\n+\n");
			for (int i = 0; i < 150; i++) sb.append((char) ('!' + 30 + random.nextInt(12)));
			sb.append('\n');
		}
		return sb.toString().getBytes(StandardCharsets.US_ASCII);
	}

	static byte [] random (int size, long seed) {
		byte [] data = new byte[size];
		new Random(seed).nextBytes(data);
		return data;
	}

	/** Text-like data with long repeats. */
	static byte [] repetitive (int size, long seed) {
		Random random = new Random(seed);
		byte [] data = new byte[size];
		for (int i = 0; i < size; i++) {
			data[i] = i > 1000 && random.nextInt(4) != 0 ? data[i - 1 - random.nextInt(1000)] : (byte) ('a' + random.nextInt(26));
		}
		return data;
	}

	// --- Gzip writers ---

	static byte [] gzip (byte [] data, int level, int strategy, int flushEvery) {
		return gzip(data, level, strategy, flushEvery, 0, null);
	}

	/**
	 * A gzip member with the given deflate settings, optionally sync flushing every
	 * flushEvery bytes (which mixes in empty stored blocks), and header flags.
	 */
	static byte [] gzip (byte [] data, int level, int strategy, int flushEvery, int flags, byte [] extra) {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		out.write(0x1f);
		out.write(0x8b);
		out.write(8);
		out.write(flags);
		out.writeBytes(new byte [] { 1, 2, 3, 4, 0, 3 });
		if ((flags & 4) != 0) {
			out.write(extra.length & 0xFF);
			out.write(extra.length >>> 8);
			out.writeBytes(extra);
		}
		if ((flags & 8) != 0) out.writeBytes("name.fastq\0".getBytes(StandardCharsets.ISO_8859_1));
		if ((flags & 16) != 0) out.writeBytes("a comment\0".getBytes(StandardCharsets.ISO_8859_1));
		if ((flags & 2) != 0) {
			CRC32 headerCrc = new CRC32();
			headerCrc.update(out.toByteArray());
			out.write((int) headerCrc.getValue());
			out.write((int) headerCrc.getValue() >>> 8);
		}

		Deflater deflater = new Deflater(level, true);
		deflater.setStrategy(strategy);
		byte [] buffer = new byte[1 << 16];
		int step = flushEvery > 0 ? flushEvery : Math.max(1, data.length);
		for (int pos = 0; pos < data.length; pos += step) {
			deflater.setInput(data, pos, Math.min(step, data.length - pos));
			int n;
			while ((n = deflater.deflate(buffer, 0, buffer.length, flushEvery > 0 ? Deflater.SYNC_FLUSH : Deflater.NO_FLUSH)) > 0) out.write(buffer, 0, n);
		}
		deflater.finish();
		while (!deflater.finished()) {
			int n = deflater.deflate(buffer);
			out.write(buffer, 0, n);
		}
		deflater.end();

		CRC32 crc = new CRC32();
		crc.update(data);
		writeInt(out, (int) crc.getValue());
		writeInt(out, data.length);
		return out.toByteArray();
	}

	static void writeInt (ByteArrayOutputStream out, int value) {
		for (int i = 0; i < 4; i++) out.write(value >>> (8 * i));
	}

	static byte [] concat (byte [] ... parts) {
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		for (byte [] p : parts) out.writeBytes(p);
		return out.toByteArray();
	}

	// --- Helpers ---

	Path write (byte [] data) throws IOException {
		Path file = Files.createTempFile(dir, "test", ".gz");
		Files.write(file, data);
		return file;
	}

	static byte [] readAll (InputStream in) throws IOException {
		try (in) {
			return in.readAllBytes();
		}
	}

	static byte [] jdk (byte [] gz) throws IOException {
		return readAll(new GZIPInputStream(new ByteArrayInputStream(gz)));
	}

	/** Checks the output equals GZIPInputStream's for various threads and chunk sizes. */
	void check (byte [] gz) throws IOException {
		byte [] expected = jdk(gz);
		Path file = write(gz);
		for (int threads : new int [] { 1, 2, 4, 8 }) {
			for (int chunk : new int [] { MIN_CHUNK, 3 * MIN_CHUNK + 12345, RapidGzipInputStream.DEFAULT_CHUNK_SIZE }) {
				byte [] actual = readAll(new RapidGzipInputStream(file, threads, chunk));
				assertArrayEquals(expected, actual, "threads=" + threads + " chunk=" + chunk);
			}
		}
	}

	// --- Tests ---

	@ParameterizedTest
	@ValueSource(ints = { 0, 1, 6, 9 })
	public void testFastqLevels (int level) throws IOException {
		check(gzip(fastq(20_000, level), level, Deflater.DEFAULT_STRATEGY, 0));
	}

	@ParameterizedTest
	@CsvSource({ "6, 1", "6, 2" })
	public void testStrategies (int level, int strategy) throws IOException {
		check(gzip(fastq(10_000, strategy), level, strategy, 0));
	}

	@Test
	public void testSyncFlushMixesBlockTypes () throws IOException {
		// Flushing often gives small fixed blocks and empty stored blocks between dynamic ones
		check(gzip(fastq(10_000, 7), 6, Deflater.DEFAULT_STRATEGY, 1000));
		check(gzip(repetitive(2_000_000, 7), 6, Deflater.DEFAULT_STRATEGY, 50_000));
	}

	@Test
	public void testChunksAreDecodedSpeculatively () throws IOException {
		byte [] data = fastq(50_000, 14);
		try (RapidGzipInputStream in = new RapidGzipInputStream(write(gzip(data, 6, Deflater.DEFAULT_STRATEGY, 0)), 4, MIN_CHUNK)) {
			assertArrayEquals(data, in.readAllBytes());
			// About 4MB compressed in 64KB chunks: all but a few should be decoded in parallel
			assertTrue(in.speculativeChunks > 40, "speculative chunks: " + in.speculativeChunks);
			assertTrue(in.sequentialDecodes < 5, "sequential decodes: " + in.sequentialDecodes);
		}
	}

	@Test
	public void testFallsBackToSequentialDecoding () throws IOException {
		// Tiny fixed and stored blocks everywhere make speculative starts miss the real boundary
		byte [] data = fastq(20_000, 15);
		try (RapidGzipInputStream in = new RapidGzipInputStream(write(gzip(data, 6, Deflater.DEFAULT_STRATEGY, 100)), 4, MIN_CHUNK)) {
			assertArrayEquals(data, in.readAllBytes());
			assertTrue(in.sequentialDecodes > 0, "sequential decodes: " + in.sequentialDecodes);
		}
	}

	@Test
	public void testIncompressibleData () throws IOException {
		check(gzip(random(1_000_000, 1), 6, Deflater.DEFAULT_STRATEGY, 0));
	}

	@Test
	public void testRepetitiveData () throws IOException {
		check(gzip(repetitive(3_000_000, 2), 9, Deflater.DEFAULT_STRATEGY, 0));
	}

	@Test
	public void testSmallInputs () throws IOException {
		check(gzip(new byte[0], 6, Deflater.DEFAULT_STRATEGY, 0));
		check(gzip(new byte [] { 'A' }, 6, Deflater.DEFAULT_STRATEGY, 0));
		check(gzip("@r\nACGT\n+\nIIII\n".getBytes(StandardCharsets.US_ASCII), 6, Deflater.DEFAULT_STRATEGY, 0));
	}

	@Test
	public void testMultipleMembers () throws IOException {
		check(concat(
				gzip(fastq(5_000, 1), 6, Deflater.DEFAULT_STRATEGY, 0),
				gzip(new byte[0], 6, Deflater.DEFAULT_STRATEGY, 0),
				gzip(fastq(5_000, 2), 1, Deflater.DEFAULT_STRATEGY, 0),
				gzip(random(200_000, 3), 6, Deflater.DEFAULT_STRATEGY, 0),
				gzip(fastq(100, 4), 9, Deflater.DEFAULT_STRATEGY, 0)));
	}

	@Test
	public void testManySmallMembers () throws IOException {
		// Like BGZF: many members of up to 64KB
		byte [] data = fastq(20_000, 5);
		ByteArrayOutputStream out = new ByteArrayOutputStream();
		for (int pos = 0; pos < data.length; pos += 65280) {
			out.writeBytes(gzip(java.util.Arrays.copyOfRange(data, pos, Math.min(data.length, pos + 65280)), 6, Deflater.DEFAULT_STRATEGY, 0, 4, new byte [] { 'B', 'C', 2, 0, 0, 0 }));
		}
		check(out.toByteArray());
	}

	@Test
	public void testHeaderFlags () throws IOException {
		byte [] data = fastq(2_000, 6);
		check(gzip(data, 6, Deflater.DEFAULT_STRATEGY, 0, 4 | 8 | 16 | 2, new byte [] { 1, 2, 3 }));
	}

	@Test
	public void testCorruptHeaderCrc () throws IOException {
		byte [] gz = gzip(fastq(100, 6), 6, Deflater.DEFAULT_STRATEGY, 0, 2, null);
		gz[10] ^= 1;
		assertThrows(ZipException.class, () -> jdk(gz));
		Path file = write(gz);
		assertThrows(ZipException.class, () -> new RapidGzipInputStream(file, 2));
	}

	@Test
	public void testTrailingGarbageIsIgnored () throws IOException {
		byte [] data = fastq(2_000, 8);
		byte [] gz = concat(gzip(data, 6, Deflater.DEFAULT_STRATEGY, 0), "not gzip".getBytes(StandardCharsets.US_ASCII));
		assertArrayEquals(data, readAll(new RapidGzipInputStream(write(gz), 4, MIN_CHUNK)));
	}

	@Test
	public void testSingleByteReads () throws IOException {
		byte [] data = fastq(3_000, 9);
		try (InputStream in = new RapidGzipInputStream(write(gzip(data, 6, Deflater.DEFAULT_STRATEGY, 0)), 4, MIN_CHUNK)) {
			for (int i = 0; i < data.length; i++) assertEquals(data[i] & 0xFF, in.read(), "byte " + i);
			assertEquals(-1, in.read());
			assertEquals(-1, in.read(new byte[10], 0, 10));
		}
	}

	@Test
	public void testCorruptCrc () throws IOException {
		byte [] gz = gzip(fastq(20_000, 10), 6, Deflater.DEFAULT_STRATEGY, 0);
		gz[gz.length - 6] ^= 1;
		Path file = write(gz);
		assertThrows(ZipException.class, () -> readAll(new RapidGzipInputStream(file, 4, MIN_CHUNK)));
	}

	@Test
	public void testCorruptSize () throws IOException {
		byte [] gz = gzip(fastq(20_000, 11), 6, Deflater.DEFAULT_STRATEGY, 0);
		gz[gz.length - 3] ^= 1;
		Path file = write(gz);
		assertThrows(ZipException.class, () -> readAll(new RapidGzipInputStream(file, 4, MIN_CHUNK)));
	}

	@Test
	public void testCorruptData () throws IOException {
		byte [] gz = gzip(fastq(20_000, 12), 6, Deflater.DEFAULT_STRATEGY, 0);
		gz[gz.length / 2] ^= 0x55;
		Path file = write(gz);
		assertThrows(IOException.class, () -> readAll(new RapidGzipInputStream(file, 4, MIN_CHUNK)));
	}

	@Test
	public void testTruncated () throws IOException {
		byte [] gz = gzip(fastq(20_000, 13), 6, Deflater.DEFAULT_STRATEGY, 0);
		for (int cut : new int [] { gz.length - 1, gz.length - 8, gz.length / 2, 20 }) {
			Path file = write(java.util.Arrays.copyOf(gz, cut));
			assertThrows(IOException.class, () -> readAll(new RapidGzipInputStream(file, 4, MIN_CHUNK)), "cut at " + cut);
		}
	}

	@Test
	public void testNotGzip () throws IOException {
		Path file = write("@r\nACGT\n+\nIIII\n".getBytes(StandardCharsets.US_ASCII));
		assertThrows(ZipException.class, () -> new RapidGzipInputStream(file, 2));
	}

	@Test
	public void testEmptyFile () throws IOException {
		Path file = write(new byte[0]);
		assertThrows(EOFException.class, () -> new RapidGzipInputStream(file, 2));
	}

}
