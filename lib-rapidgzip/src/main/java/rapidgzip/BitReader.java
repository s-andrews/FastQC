package rapidgzip;

import java.io.EOFException;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;

/**
 * Least significant bit first reader over compressed data, as deflate needs.
 * Data is loaded on demand from a {@link Source} by absolute file offset, so a
 * reader can start at any bit of the file.
 *
 * Bits are buffered in a 64 bit word, refilled 8 bytes at a time with one
 * unaligned little endian load. Past the end of the data it supplies zero
 * bytes instead of failing in the middle of a symbol; callers check
 * {@link #pastEnd()} at block and trailer boundaries. Not thread safe: one
 * reader per decode.
 */
final class BitReader {

	/** Random access to the compressed bytes. */
	interface Source {
		/** Reads up to {@code length} bytes at {@code offset}, fewer only at the end of the data. */
		byte [] read (long offset, int length) throws IOException;
	}

	private static final VarHandle LONG_LE = MethodHandles.byteArrayViewVarHandle(long[].class, ByteOrder.LITTLE_ENDIAN);

	// Reading this many zero bytes beyond the end of the data means a truncated stream
	private static final int MAX_PADDING = 64;

	private final Source source;
	private final int segmentSize;

	private byte [] data = new byte[0];
	// File offset of data[0]
	private long base;
	private int limit;
	private int pos;
	private boolean endOfData;
	// Zero bytes supplied past the end of the data
	private int padding;

	// Bit buffer, the next bit to read is the lowest. Package private so the decoder's
	// hot loops can read it directly after a refill.
	long bits;
	int count;

	BitReader (Source source, int segmentSize) {
		this.source = source;
		this.segmentSize = segmentSize;
	}

	/** Absolute position of the next bit to read. */
	long position () {
		return (base + pos + padding) * 8 - count;
	}

	/** True once the reader has gone past the last byte of the data. */
	boolean pastEnd () {
		return padding > 0 && position() > (base + limit) * 8;
	}

	void seek (long bitPosition) throws IOException {
		long byteOffset = bitPosition >>> 3;
		bits = 0;
		count = 0;
		padding = 0;
		if (byteOffset >= base && byteOffset <= base + limit) {
			pos = (int) (byteOffset - base);
		}
		else {
			base = byteOffset;
			data = source.read(byteOffset, segmentSize);
			limit = data.length;
			pos = 0;
			endOfData = limit < segmentSize;
		}
		int skip = (int) (bitPosition & 7);
		if (skip > 0) {
			refill();
			consume(skip);
		}
	}

	/**
	 * Makes at least 56 bits available. The fast path assumes count < 64, which
	 * holds as every caller consumes bits between refills (the slow path can
	 * leave count at 64).
	 */
	void refill () throws IOException {
		if (pos + 8 <= limit) {
			long word = (long) LONG_LE.get(data, pos);
			bits |= word << count;
			int bytes = (63 - count) >>> 3;
			pos += bytes;
			count += bytes << 3;
			bits &= (1L << count) - 1;
		}
		else {
			refillSlow();
		}
	}

	private void refillSlow () throws IOException {
		while (count <= 56) {
			if (pos == limit && !loadMore()) {
				// Past the end: supply zeros, which decoding fails on soon enough
				if (++padding > MAX_PADDING) throw new EOFException("Unexpected end of compressed data");
				count += 8;
				continue;
			}
			bits |= (data[pos++] & 0xFFL) << count;
			count += 8;
		}
	}

	// Appends the next segment to the unread bytes (earlier bytes are dropped)
	private boolean loadMore () throws IOException {
		if (endOfData) return false;
		byte [] more = source.read(base + limit, segmentSize);
		if (more.length < segmentSize) endOfData = true;
		if (more.length == 0) return false;
		int keep = limit - pos;
		byte [] joined = new byte[keep + more.length];
		System.arraycopy(data, pos, joined, 0, keep);
		System.arraycopy(more, 0, joined, keep, more.length);
		base += pos;
		data = joined;
		limit = joined.length;
		pos = 0;
		return true;
	}

	/** Lowest n bits, needs n <= count. */
	int peek (int n) {
		return (int) (bits & ((1L << n) - 1));
	}

	void consume (int n) {
		bits >>>= n;
		count -= n;
	}

	int read (int n) throws IOException {
		if (count < n) refill();
		int value = peek(n);
		consume(n);
		return value;
	}

	void alignToByte () {
		consume(count & 7);
	}

	/** Reads a byte, the reader must be byte aligned. */
	int readByte () throws IOException {
		return read(8);
	}

	/** Copies n bytes into dest, the reader must be byte aligned. */
	void readBytes (byte [] dest, int offset, int n) throws IOException {
		// Drain whole bytes left in the bit buffer first
		while (n > 0 && count >= 8) {
			dest[offset++] = (byte) peek(8);
			consume(8);
			n--;
		}
		while (n > 0) {
			if (pos == limit && !loadMore()) throw new EOFException("Unexpected end of compressed data");
			int chunk = Math.min(n, limit - pos);
			System.arraycopy(data, pos, dest, offset, chunk);
			pos += chunk;
			offset += chunk;
			n -= chunk;
		}
	}

	/** True if there's no byte left after the current (byte aligned) position. */
	boolean atEnd () throws IOException {
		if (count - padding * 8 >= 8) return false;
		if (padding > 0) return true;
		if (pos < limit) return false;
		return !loadMore();
	}

}
