package rapidgzip;

import java.io.EOFException;
import java.io.IOException;
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
 * {@link Chunk}. Output before {@link #lastMarker} may contain markers, output
 * after it doesn't. In FASTQ about 55% of the output refers back into the
 * window, so most chunks never leave marker mode.
 *
 * Stops only at block boundaries (checked before each block, after the first),
 * so the end bit of a chunk is always a real boundary when its start was one.
 * Member boundaries are decoded inline: trailer, then the next gzip header.
 * One instance per decode, used by a single thread.
 */
final class DeflateDecoder {

	static final int WINDOW = 32768;

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

	// Marker mode output, starting with the WINDOW markers for the unknown window
	private char [] chars;
	private int charPos;
	// Output before this position may contain markers
	private int lastMarker;
	private boolean markerMode;

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

	DeflateDecoder (BitReader in) {
		this.in = in;
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
			chars = new char[this.sizeHint + WINDOW];
			for (int i = 0; i < WINDOW; i++) chars[i] = (char) (256 + i);
			charPos = WINDOW;
			lastMarker = WINDOW;
		}
		else {
			markerMode = false;
			bytes = new byte[this.sizeHint + window.length];
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
	 * Byte mode is cheaper: half the memory and no marker tracking per copy.
	 */
	private void switchToBytes () {
		markerMode = false;
		bytes = new byte[Math.max(sizeHint - (charPos - WINDOW), 1 << 16) + WINDOW];
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
			if (charPos - lastMarker >= WINDOW) switchToBytes();
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
		if (charPos + n > chars.length) chars = Arrays.copyOf(chars, Math.max(chars.length + (chars.length >> 1), charPos + n));
	}

	private void ensureBytes (int n) {
		if (bytePos + n > bytes.length) bytes = Arrays.copyOf(bytes, Math.max(bytes.length + (bytes.length >> 1), bytePos + n));
	}

	/**
	 * Marker mode decoding. Returns true at the end of the block, or false
	 * after switching to byte mode to continue the block there.
	 *
	 * Hot loop: state is kept in locals and written back in the finally block.
	 * A match records where it ended if it copied any marker (OR of the values
	 * copied >= 256), which is all lastMarker needs; literals never are markers.
	 * Each refill gives >= 56 bits, enough for a literal/length code, its extra
	 * bits, a distance code and its extra bits (at most 48).
	 */
	private boolean decodeChars (int [] lit, int [] dist) throws IOException {
		char [] out = chars;
		int op = charPos;
		int marker = lastMarker;
		try {
			while (true) {
				if (op - marker >= WINDOW) {
					charPos = op;
					lastMarker = marker;
					chars = out;
					switchToBytes();
					return false;
				}
				if (op + MAX_MATCH > out.length) {
					charPos = op;
					ensureChars(MAX_MATCH);
					out = chars;
				}
				in.refill();
				long b = in.bits;
				int e = lit[(int) b & LIT_MASK];
				if (e < 0) e = lit[((e >>> 4) & 0x7FFFFFF) + ((int) (b >>> LIT_BITS) & ((1 << (e & 15)) - 1))];
				int n = e & 15;
				b >>>= n;
				int used = n;
				int symbol = e >>> 4;
				if (symbol < 256) {
					out[op++] = (char) symbol;
					in.consume(used);
					continue;
				}
				if (symbol == 256) {
					in.consume(used);
					return true;
				}
				symbol -= 257;
				if (symbol >= 29) throw new DataFormatError("Invalid literal/length code");
				int extra = LENGTH_EXTRA[symbol];
				int length = LENGTH_BASE[symbol] + ((int) b & ((1 << extra) - 1));
				b >>>= extra;
				used += extra;

				e = dist[(int) b & DIST_MASK];
				if (e < 0) e = dist[((e >>> 4) & 0x7FFFFFF) + ((int) (b >>> DIST_BITS) & ((1 << (e & 15)) - 1))];
				n = e & 15;
				b >>>= n;
				used += n;
				symbol = e >>> 4;
				if (symbol >= 30) throw new DataFormatError("Invalid distance code");
				extra = DIST_EXTRA[symbol];
				int distance = DIST_BASE[symbol] + ((int) b & ((1 << extra) - 1));
				used += extra;
				in.consume(used);

				// The WINDOW markers at the start of the output make any distance <= 32KB valid
				int from = op - distance;
				if (from < 0) throw new DataFormatError("Invalid distance too far back");
				// Element by element, as overlapping copies (distance < length) must repeat the pattern
				int seen = 0;
				for (int k = 0; k < length; k++) {
					char c = out[from + k];
					out[op + k] = c;
					seen |= c;
				}
				op += length;
				if (seen >= 256) marker = op;
			}
		}
		finally {
			if (markerMode) {
				charPos = op;
				lastMarker = marker;
				chars = out;
			}
		}
	}

	/** Byte mode decoding to the end of the block: the window is known, as in a plain inflater. */
	private void decodeBytes (int [] lit, int [] dist) throws IOException {
		byte [] out = bytes;
		int op = bytePos;
		try {
			while (true) {
				if (op + MAX_MATCH > out.length) {
					bytePos = op;
					ensureBytes(MAX_MATCH);
					out = bytes;
				}
				in.refill();
				long b = in.bits;
				int e = lit[(int) b & LIT_MASK];
				if (e < 0) e = lit[((e >>> 4) & 0x7FFFFFF) + ((int) (b >>> LIT_BITS) & ((1 << (e & 15)) - 1))];
				int n = e & 15;
				b >>>= n;
				int used = n;
				int symbol = e >>> 4;
				if (symbol < 256) {
					out[op++] = (byte) symbol;
					in.consume(used);
					continue;
				}
				if (symbol == 256) {
					in.consume(used);
					return;
				}
				symbol -= 257;
				if (symbol >= 29) throw new DataFormatError("Invalid literal/length code");
				int extra = LENGTH_EXTRA[symbol];
				int length = LENGTH_BASE[symbol] + ((int) b & ((1 << extra) - 1));
				b >>>= extra;
				used += extra;

				e = dist[(int) b & DIST_MASK];
				if (e < 0) e = dist[((e >>> 4) & 0x7FFFFFF) + ((int) (b >>> DIST_BITS) & ((1 << (e & 15)) - 1))];
				n = e & 15;
				b >>>= n;
				used += n;
				symbol = e >>> 4;
				if (symbol >= 30) throw new DataFormatError("Invalid distance code");
				extra = DIST_EXTRA[symbol];
				int distance = DIST_BASE[symbol] + ((int) b & ((1 << extra) - 1));
				used += extra;
				in.consume(used);

				// out starts with the (up to 32KB) window, so this rejects what zlib rejects
				int from = op - distance;
				if (from < 0) throw new DataFormatError("Invalid distance too far back");
				// arraycopy only when source and destination don't overlap
				if (distance >= length) {
					System.arraycopy(out, from, out, op, length);
				}
				else {
					for (int k = 0; k < length; k++) out[op + k] = out[from + k];
				}
				op += length;
			}
		}
		finally {
			bytePos = op;
			bytes = out;
		}
	}

}
