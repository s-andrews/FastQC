package rapidgzip;

import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;

/**
 * Finds likely starts of non-final dynamic Huffman deflate blocks by checking
 * every bit offset for a valid block header: block type, symbol counts, a
 * complete precode, valid code length runs, and valid literal/length and
 * distance codes including an end-of-block code. Most compressors (zlib, igzip,
 * libdeflate...) emit almost only such blocks.
 */
final class BlockFinder {

	private static final VarHandle LONG_LE = MethodHandles.byteArrayViewVarHandle(long[].class, ByteOrder.LITTLE_ENDIAN);

	// Indexed by the first 13 bits at an offset: BFINAL=0, BTYPE=2, HLIT <= 29, HDIST <= 29
	private static final boolean [] HEADER_OK = new boolean[1 << 13];
	static {
		for (int h = 0; h < HEADER_OK.length; h++) {
			HEADER_OK[h] = (h & 7) == 0b100 && ((h >>> 3) & 31) <= 29 && ((h >>> 8) & 31) <= 29;
		}
	}

	private final byte [] precodeLengths = new byte[19];
	private final byte [] codeLengths = new byte[320];
	private final int [] precodeCounts = new int[8];

	/**
	 * Returns the first bit offset in [fromBit, toBit) of data which looks like a
	 * dynamic block start, or -1. Bit offsets are relative to data[0], and only
	 * offsets with at least maxHeaderBytes of data after them are checked.
	 */
	long find (byte [] data, int length, long fromBit, long toBit) {
		// A dynamic block header is at most 3+14+57+(286+30)*7 bits
		long lastBit = Math.min(toBit, (long) (length - 8) * 8);
		for (long bit = fromBit; bit < lastBit; bit++) {
			int bytePos = (int) (bit >>> 3);
			long word = (long) LONG_LE.get(data, bytePos) >>> (bit & 7);
			if (!HEADER_OK[(int) word & 0x1FFF]) continue;
			if (isDynamicHeader(data, length, bit)) return bit;
		}
		return -1;
	}

	private boolean isDynamicHeader (byte [] data, int length, long bit) {
		long pos = bit + 3;
		long limit = (long) length * 8;
		if (pos + 14 + 19*3 > limit) return false;
		int literals = bits(data, pos, 5) + 257;
		int distances = bits(data, pos + 5, 5) + 1;
		int precodes = bits(data, pos + 10, 4) + 4;
		pos += 14;

		// The precode must be complete
		java.util.Arrays.fill(precodeLengths, (byte) 0);
		java.util.Arrays.fill(precodeCounts, 0);
		for (int i = 0; i < precodes; i++) {
			int l = bits(data, pos, 3);
			pos += 3;
			precodeLengths[DeflateDecoder.PRECODE_ORDER[i]] = (byte) l;
			precodeCounts[l]++;
		}
		int left = 1;
		for (int l = 1; l <= 7; l++) {
			left = (left << 1) - precodeCounts[l];
			if (left < 0) return false;
		}
		if (left != 0) return false;

		int [] precode;
		try {
			precode = Huffman.build(precodeLengths, 0, 19, 7, 19, false);
		}
		catch (DataFormatError e) {
			return false;
		}

		int total = literals + distances;
		int i = 0;
		while (i < total) {
			if (pos + 7 + 7 > limit) return false;
			int entry = precode[bits(data, pos, 7)];
			pos += entry & 15;
			int symbol = entry >>> 4;
			if (symbol < 16) {
				codeLengths[i++] = (byte) symbol;
				continue;
			}
			int repeat;
			byte value = 0;
			if (symbol == 16) {
				if (i == 0) return false;
				value = codeLengths[i-1];
				repeat = 3 + bits(data, pos, 2);
				pos += 2;
			}
			else if (symbol == 17) {
				repeat = 3 + bits(data, pos, 3);
				pos += 3;
			}
			else {
				repeat = 11 + bits(data, pos, 7);
				pos += 7;
			}
			if (i + repeat > total) return false;
			java.util.Arrays.fill(codeLengths, i, i + repeat, value);
			i += repeat;
		}
		if (codeLengths[256] == 0) return false;
		return Huffman.isValid(codeLengths, 0, literals, false) && Huffman.isValid(codeLengths, literals, distances, true);
	}

	private static int bits (byte [] data, long bit, int n) {
		int bytePos = (int) (bit >>> 3);
		int shift = (int) (bit & 7);
		int value = data[bytePos] & 0xFF;
		if (shift + n > 8) value |= (data[bytePos+1] & 0xFF) << 8;
		if (shift + n > 16) value |= (data[bytePos+2] & 0xFF) << 16;
		return (value >>> shift) & ((1 << n) - 1);
	}

}
