package rapidgzip;

import java.io.EOFException;
import java.io.IOException;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;
import java.util.Arrays;

/**
 * Decodes deflate blocks, and the gzip framing between members, from a block
 * boundary up to the first block boundary at or after a stop position.
 *
 * Without the preceding window (when decoding speculatively from a guessed
 * boundary) output starts as 16 bit values with markers for the unknown window
 * bytes, and switches to plain bytes once the last 32KB of output contain no
 * marker, as nothing before them can be referenced any more.
 *
 * Marker mode: the output buffer starts with WINDOW values 256..256+32767
 * standing for the unknown window, so a back reference of any valid distance
 * just copies values, markers included, and is resolved later by
 * {@link Chunk}. Whether the last WINDOW values still hold a marker is checked
 * every MARKER_CHECK values (see markersInWindow); once none do, output switches
 * to bytes. In FASTQ about 55% of the output refers back into the window, so
 * most chunks never leave marker mode.
 *
 * Stops only at block boundaries (checked before each block, after the first),
 * so the end bit of a chunk is always a real boundary when its start was one.
 * Member boundaries are decoded inline: trailer, then the next gzip header.
 * One instance per decode, used by a single thread.
 */
final class DeflateDecoder {

	static final int WINDOW = 32768;

	private static final VarHandle LONG_LE = MethodHandles.byteArrayViewVarHandle(long[].class, ByteOrder.LITTLE_ENDIAN);
	// Bits needed to decode a length/distance pair: 15 + 5 + 15 + 13. The decoding
	// loops refill only below this, to 56-63 bits, so a whole symbol decodes without
	// further checks and several literals share one refill.
	private static final int MAX_SYMBOL_BITS = 48;
	// Below this, back references are copied with a loop rather than System.arraycopy:
	// for the short matches common in FASTQ the call costs more than the copy
	private static final int SHORT_COPY = 32;

	// Primary table sizes (see Huffman); longer codes go through a second level.
	// Out of range symbols for invalid bit patterns, rejected when decoded.
	private static final int LIT_BITS = 10;
	private static final int LIT_MASK = (1 << LIT_BITS) - 1;
	private static final int DIST_BITS = 8;
	private static final int DIST_MASK = (1 << DIST_BITS) - 1;
	private static final int INVALID_LIT = 300;
	private static final int INVALID_DIST = 31;
	private static final int INVALID_PRECODE = 19;
	private static final int MAX_MATCH = 258;

	private static final int [] LENGTH_BASE = { 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19, 23, 27, 31, 35, 43, 51, 59, 67, 83, 99, 115, 131, 163, 195, 227, 258 };
	private static final int [] LENGTH_EXTRA = { 0, 0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3, 4, 4, 4, 4, 5, 5, 5, 5, 0 };
	private static final int [] DIST_BASE = { 1, 2, 3, 4, 5, 7, 9, 13, 17, 25, 33, 49, 65, 97, 129, 193, 257, 385, 513, 769, 1025, 1537, 2049, 3073, 4097, 6145, 8193, 12289, 16385, 24577 };
	private static final int [] DIST_EXTRA = { 0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 5, 5, 6, 6, 7, 7, 8, 8, 9, 9, 10, 10, 11, 11, 12, 12, 13, 13 };
	static final int [] PRECODE_ORDER = { 16, 17, 18, 0, 8, 7, 9, 6, 10, 5, 11, 4, 12, 3, 13, 2, 14, 1, 15 };

	private static final int [] FIXED_LIT;
	private static final int [] FIXED_DIST;
	static {
		byte [] lengths = new byte[320];
		Arrays.fill(lengths, 0, 144, (byte) 8);
		Arrays.fill(lengths, 144, 256, (byte) 9);
		Arrays.fill(lengths, 256, 280, (byte) 7);
		Arrays.fill(lengths, 280, 288, (byte) 8);
		Arrays.fill(lengths, 288, 320, (byte) 5);
		try {
			FIXED_LIT = Huffman.build(lengths, 0, 288, LIT_BITS, INVALID_LIT, false);
			FIXED_DIST = Huffman.build(lengths, 288, 32, DIST_BITS, INVALID_DIST, false);
		}
		catch (DataFormatError e) {
			throw new IllegalStateException(e);
		}
	}

	private final BitReader in;
	// Output buffers come from, and outgrown ones go back to, the stream's pool
	private final Buffers buffers;

	// Marker mode output, starting with the WINDOW markers for the unknown window
	private char [] chars;
	private int charPos;
	private boolean markerMode;
	// Output position at which to next check whether markers are gone (see markersInWindow)
	private int nextMarkerCheck;

	// Byte mode output, starting with the known window
	private byte [] bytes;
	private int bytePos;
	private int bytesFrom;

	private long [] memberEnds = new long[0];
	private int [] memberCrcs = new int[0];
	private int [] memberSizes = new int[0];

	// Expected output size, to allocate the output buffer once (growing copies it)
	private int sizeHint;

	private final byte [] codeLengths = new byte[320];
	private final byte [] precodeLengths = new byte[19];

	DeflateDecoder (BitReader in, Buffers buffers) {
		this.in = in;
		this.buffers = buffers;
	}

	/**
	 * Decodes from the block boundary at startBit to the first block boundary at
	 * or after stopBit, or the end of the data.
	 *
	 * @param window the up to 32KB of output preceding startBit, or null if unknown
	 * @throws DataFormatError on invalid data, which for a guessed startBit means it wasn't a block start
	 * @throws EOFException if the data ends before a block boundary at or after stopBit
	 */
	Chunk decode (long startBit, long stopBit, byte [] window, int sizeHint) throws IOException {
		this.sizeHint = Math.max(sizeHint, 1 << 16);
		if (window == null) {
			// Markers often last for the whole chunk (FASTQ), so size for it
			markerMode = true;
			chars = buffers.chars(this.sizeHint + WINDOW);
			for (int i = 0; i < WINDOW; i++) chars[i] = (char) (256 + i);
			charPos = WINDOW;
			// Until WINDOW values are output, the window still holds initial markers
			nextMarkerCheck = 2 * WINDOW;
		}
		else {
			markerMode = false;
			bytes = buffers.bytes(this.sizeHint + window.length);
			System.arraycopy(window, 0, bytes, 0, window.length);
			bytePos = bytesFrom = window.length;
		}

		in.seek(startBit);
		boolean endOfData = false;
		boolean first = true;
		while (true) {
			// At least one block, so a sequential decode always makes progress
			if (!first && in.position() >= stopBit) break;
			first = false;

			int header = in.read(3);
			switch (header >>> 1) {
				case 0: storedBlock(); break;
				case 1: huffmanBlock(FIXED_LIT, FIXED_DIST); break;
				case 2: dynamicBlock(); break;
				default: throw new DataFormatError("Invalid block type");
			}
			// BitReader supplies zeros past the end of the data: a block that only
			// decoded thanks to them is a truncated stream
			if (in.pastEnd()) throw new EOFException("Unexpected end of compressed data");

			if ((header & 1) != 0) {
				// Last block of a gzip member: the trailer, then maybe another member
				in.alignToByte();
				int crc = readInt32();
				int size = readInt32();
				if (in.pastEnd()) throw new EOFException("Unexpected end of compressed data");
				// Checked by the reader, which sees the output in order. The output decoded so
				// far stays the window of the next member, which valid data never references.
				addMemberEnd(crc, size);
				if (in.atEnd() || !GzipHeader.skip(in)) {
					// Like GZIPInputStream, anything after the last member that isn't one is ignored
					endOfData = true;
					break;
				}
			}
		}

		long endBit = in.position();
		// Output is chars[WINDOW, charPos) if any marker mode output, then bytes[bytesFrom, bytePos)
		if (markerMode) {
			return new Chunk(startBit, endBit, endOfData, chars, WINDOW, charPos, new byte[0], 0, 0, memberEnds, memberCrcs, memberSizes);
		}
		int markersTo = chars == null ? WINDOW : charPos;
		return new Chunk(startBit, endBit, endOfData, chars, WINDOW, markersTo, bytes, bytesFrom, bytePos, memberEnds, memberCrcs, memberSizes);
	}

	private int readInt32 () throws IOException {
		return in.readByte() | in.readByte() << 8 | in.readByte() << 16 | in.readByte() << 24;
	}

	private long outputLength () {
		long markerOutput = chars == null ? 0 : charPos - WINDOW;
		long byteOutput = markerMode ? 0 : bytePos - bytesFrom;
		return markerOutput + byteOutput;
	}

	private void addMemberEnd (int crc, int size) {
		int n = memberEnds.length;
		memberEnds = Arrays.copyOf(memberEnds, n+1);
		memberCrcs = Arrays.copyOf(memberCrcs, n+1);
		memberSizes = Arrays.copyOf(memberSizes, n+1);
		memberEnds[n] = outputLength();
		memberCrcs[n] = crc;
		memberSizes[n] = size;
	}

	/**
	 * Switches to byte mode, the last WINDOW values being free of markers. They
	 * are copied as the byte mode window, and stay in chars as output too, so the
	 * chunk's output is chars[WINDOW, charPos) then bytes[WINDOW, bytePos).
	 * Byte mode is cheaper: half the memory and no periodic marker scans.
	 */
	private void switchToBytes () {
		markerMode = false;
		bytes = buffers.bytes(Math.max(sizeHint - (charPos - WINDOW), 1 << 16) + WINDOW);
		for (int i = 0; i < WINDOW; i++) bytes[i] = (byte) chars[charPos - WINDOW + i];
		bytePos = bytesFrom = WINDOW;
	}

	private void storedBlock () throws IOException {
		in.alignToByte();
		int length = in.read(16);
		int complement = in.read(16);
		if (length != (~complement & 0xFFFF)) throw new DataFormatError("Invalid stored block lengths");
		if (markerMode) {
			byte [] data = new byte[length];
			in.readBytes(data, 0, length);
			ensureChars(length);
			for (int i = 0; i < length; i++) chars[charPos++] = (char) (data[i] & 0xFF);
			// Stored data is literal, so it can push the markers out of the window
			if (charPos >= nextMarkerCheck && !markersInWindow()) switchToBytes();
		}
		else {
			ensureBytes(length);
			in.readBytes(bytes, bytePos, length);
			bytePos += length;
		}
	}

	private void dynamicBlock () throws IOException {
		in.refill();
		int literals = in.peek(5) + 257;
		in.consume(5);
		int distances = in.peek(5) + 1;
		in.consume(5);
		int precodes = in.peek(4) + 4;
		in.consume(4);
		if (literals > 286 || distances > 30) throw new DataFormatError("Too many length or distance symbols");

		Arrays.fill(precodeLengths, (byte) 0);
		for (int i = 0; i < precodes; i++) {
			precodeLengths[PRECODE_ORDER[i]] = (byte) in.read(3);
		}
		// Only complete precodes are valid, so all codes fit in the 7 bit table.
		// Same validity rules as zlib, so the same streams are rejected.
		int [] precode = Huffman.build(precodeLengths, 0, 19, 7, INVALID_PRECODE, false);
		if (!Huffman.isValid(precodeLengths, 0, 19, false) || hasSingleCode(precodeLengths)) throw new DataFormatError("Invalid code lengths set");

		int total = literals + distances;
		int i = 0;
		while (i < total) {
			in.refill();
			int entry = precode[in.peek(7)];
			in.consume(entry & 15);
			int symbol = entry >>> 4;
			if (symbol < 16) {
				codeLengths[i++] = (byte) symbol;
				continue;
			}
			int repeat;
			byte value = 0;
			if (symbol == 16) {
				if (i == 0) throw new DataFormatError("Invalid bit length repeat");
				value = codeLengths[i-1];
				repeat = 3 + in.peek(2);
				in.consume(2);
			}
			else if (symbol == 17) {
				repeat = 3 + in.peek(3);
				in.consume(3);
			}
			else if (symbol == 18) {
				repeat = 11 + in.peek(7);
				in.consume(7);
			}
			else {
				throw new DataFormatError("Invalid code lengths set");
			}
			if (i + repeat > total) throw new DataFormatError("Invalid bit length repeat");
			Arrays.fill(codeLengths, i, i + repeat, value);
			i += repeat;
		}
		if (codeLengths[256] == 0) throw new DataFormatError("Invalid code -- missing end-of-block");

		int [] lit = Huffman.build(codeLengths, 0, literals, LIT_BITS, INVALID_LIT, false);
		int [] dist = Huffman.build(codeLengths, literals, distances, DIST_BITS, INVALID_DIST, true);
		huffmanBlock(lit, dist);
	}

	private static boolean hasSingleCode (byte [] lengths) {
		int codes = 0;
		for (byte l : lengths) if (l != 0) codes++;
		return codes == 1;
	}

	private void huffmanBlock (int [] lit, int [] dist) throws IOException {
		if (markerMode && decodeChars(lit, dist)) return;
		decodeBytes(lit, dist);
	}

	private void ensureChars (int n) {
		if (charPos + n > chars.length) {
			char [] old = chars;
			chars = Arrays.copyOf(old, Math.max(old.length + (old.length >> 1), charPos + n));
			buffers.release(old);
		}
	}

	private void ensureBytes (int n) {
		if (bytePos + n > bytes.length) {
			byte [] old = bytes;
			bytes = Arrays.copyOf(old, Math.max(old.length + (old.length >> 1), bytePos + n));
			buffers.release(old);
		}
	}

	/**
	 * True if the last WINDOW values of output contain a marker. Checked at intervals
	 * of output rather than tracking markers through copies, which kept back references
	 * from being plain array copies. The interval doubles at each check, from 128K up
	 * to MAX_MARKER_CHECK values, as in FASTQ markers mostly last the whole chunk.
	 * Switching to bytes later than possible leaves the output unchanged (the same
	 * values, resolved later), only holds more output as 16 bit values.
	 */
	private boolean markersInWindow () {
		char [] out = chars;
		int seen = 0;
		for (int i = charPos - WINDOW; i < charPos; i++) seen |= out[i];
		// Markers often last for the whole chunk, so check less and less often
		markerCheckInterval = Math.min(markerCheckInterval * 2, MAX_MARKER_CHECK);
		nextMarkerCheck = charPos + markerCheckInterval;
		return seen >= 256;
	}

	private static final int MAX_MARKER_CHECK = 1 << 20;
	private int markerCheckInterval = 1 << 16;

	/**
	 * Marker mode decoding. Returns true at the end of the block, or false
	 * after switching to byte mode to continue the block there. The WINDOW
	 * markers at the start of the output make any distance up to 32KB valid.
	 *
	 * The bit buffer is kept in locals (see BitReader) rather than refilled and
	 * consumed through the reader for each symbol, and written back on the way out.
	 * With the short copies, marker check backoff and pooled buffers, decompression
	 * alone went from 1.53 s to 1.31 s on 4 threads (510 MB AVITI), FastQC CPU -6%
	 * (AVITI) / -2% (WES).
	 */
	private boolean decodeChars (int [] lit, int [] dist) throws IOException {
		char [] out = chars;
		int op = charPos;
		BitReader in = this.in;
		long bits = in.bits;
		int count = in.count;
		byte [] data = in.data;
		int pos = in.pos;
		int limit = in.limit;
		try {
			while (true) {
				if (op >= nextMarkerCheck) {
					charPos = op;
					if (!markersInWindow()) {
						in.bits = bits;
						in.count = count;
						in.pos = pos;
						switchToBytes();
						return false;
					}
				}
				if (op + MAX_MATCH > out.length) {
					charPos = op;
					ensureChars(MAX_MATCH);
					out = chars;
				}
				if (count < MAX_SYMBOL_BITS) {
					if (pos + 8 <= limit) {
						bits |= (long) LONG_LE.get(data, pos) << count;
						int bytes = (63 - count) >>> 3;
						pos += bytes;
						count += bytes << 3;
						bits &= (1L << count) - 1;
					}
					else {
						in.bits = bits;
						in.count = count;
						in.pos = pos;
						in.refill();
						bits = in.bits;
						count = in.count;
						data = in.data;
						pos = in.pos;
						limit = in.limit;
					}
				}
				int e = lit[(int) bits & LIT_MASK];
				if (e < 0) e = lit[((e >>> 4) & 0x7FFFFFF) + ((int) (bits >>> LIT_BITS) & ((1 << (e & 15)) - 1))];
				int n = e & 15;
				bits >>>= n;
				count -= n;
				int symbol = e >>> 4;
				if (symbol < 256) {
					out[op++] = (char) symbol;
					continue;
				}
				if (symbol == 256) {
					return true;
				}
				symbol -= 257;
				if (symbol >= 29) throw new DataFormatError("Invalid literal/length code");
				int extra = LENGTH_EXTRA[symbol];
				int length = LENGTH_BASE[symbol] + ((int) bits & ((1 << extra) - 1));
				bits >>>= extra;
				count -= extra;

				e = dist[(int) bits & DIST_MASK];
				if (e < 0) e = dist[((e >>> 4) & 0x7FFFFFF) + ((int) (bits >>> DIST_BITS) & ((1 << (e & 15)) - 1))];
				n = e & 15;
				bits >>>= n;
				count -= n;
				symbol = e >>> 4;
				if (symbol >= 30) throw new DataFormatError("Invalid distance code");
				extra = DIST_EXTRA[symbol];
				int distance = DIST_BASE[symbol] + ((int) bits & ((1 << extra) - 1));
				bits >>>= extra;
				count -= extra;

				int from = op - distance;
				if (from < 0) throw new DataFormatError("Invalid distance too far back");
				// Long copies whose source and destination don't overlap: arraycopy (an intrinsic)
				// is faster than the loop. Overlapping copies repeat the last distance values, which
				// needs the forward loop (arraycopy would copy as if through a temporary buffer).
				if (distance >= length && length > SHORT_COPY) {
					System.arraycopy(out, from, out, op, length);
				}
				else {
					for (int k = 0; k < length; k++) out[op + k] = out[from + k];
				}
				op += length;
			}
		}
		finally {
			in.bits = bits;
			in.count = count;
			in.pos = pos;
			if (markerMode) {
				charPos = op;
				chars = out;
			}
		}
	}

	/**
	 * Byte mode decoding to the end of the block, as decodeChars: the window is known,
	 * as in a plain inflater. Always returns true, at the end of the block.
	 */
	private boolean decodeBytes (int [] lit, int [] dist) throws IOException {
		byte [] out = bytes;
		int op = bytePos;
		BitReader in = this.in;
		long bits = in.bits;
		int count = in.count;
		byte [] data = in.data;
		int pos = in.pos;
		int limit = in.limit;
		try {
			while (true) {
				if (op + MAX_MATCH > out.length) {
					bytePos = op;
					ensureBytes(MAX_MATCH);
					out = bytes;
				}
				if (count < MAX_SYMBOL_BITS) {
					if (pos + 8 <= limit) {
						bits |= (long) LONG_LE.get(data, pos) << count;
						int bytes = (63 - count) >>> 3;
						pos += bytes;
						count += bytes << 3;
						bits &= (1L << count) - 1;
					}
					else {
						in.bits = bits;
						in.count = count;
						in.pos = pos;
						in.refill();
						bits = in.bits;
						count = in.count;
						data = in.data;
						pos = in.pos;
						limit = in.limit;
					}
				}
				int e = lit[(int) bits & LIT_MASK];
				if (e < 0) e = lit[((e >>> 4) & 0x7FFFFFF) + ((int) (bits >>> LIT_BITS) & ((1 << (e & 15)) - 1))];
				int n = e & 15;
				bits >>>= n;
				count -= n;
				int symbol = e >>> 4;
				if (symbol < 256) {
					out[op++] = (byte) symbol;
					continue;
				}
				if (symbol == 256) {
					return true;
				}
				symbol -= 257;
				if (symbol >= 29) throw new DataFormatError("Invalid literal/length code");
				int extra = LENGTH_EXTRA[symbol];
				int length = LENGTH_BASE[symbol] + ((int) bits & ((1 << extra) - 1));
				bits >>>= extra;
				count -= extra;

				e = dist[(int) bits & DIST_MASK];
				if (e < 0) e = dist[((e >>> 4) & 0x7FFFFFF) + ((int) (bits >>> DIST_BITS) & ((1 << (e & 15)) - 1))];
				n = e & 15;
				bits >>>= n;
				count -= n;
				symbol = e >>> 4;
				if (symbol >= 30) throw new DataFormatError("Invalid distance code");
				extra = DIST_EXTRA[symbol];
				int distance = DIST_BASE[symbol] + ((int) bits & ((1 << extra) - 1));
				bits >>>= extra;
				count -= extra;

				// out starts with the (up to 32KB) window, so this rejects what zlib rejects
				int from = op - distance;
				if (from < 0) throw new DataFormatError("Invalid distance too far back");
				// arraycopy only for long copies whose source and destination don't overlap
				if (distance >= length && length > SHORT_COPY) {
					System.arraycopy(out, from, out, op, length);
				}
				else {
					for (int k = 0; k < length; k++) out[op + k] = out[from + k];
				}
				op += length;
			}
		}
		finally {
			in.bits = bits;
			in.count = count;
			in.pos = pos;
			bytePos = op;
			bytes = out;
		}
	}

}
