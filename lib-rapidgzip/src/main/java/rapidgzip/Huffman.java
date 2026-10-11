package rapidgzip;

/**
 * Canonical Huffman decoding tables for deflate.
 *
 * A table is indexed by the next {@code primaryBits} bits of input. An entry is
 * either {@code symbol << 4 | codeLength}, or for codes longer than the primary
 * bits a negative pointer {@code 1 << 31 | subtableOffset << 4 | subtableBits}
 * to a second level indexed by the following bits. Bit patterns which aren't a
 * code (only possible in the incomplete codes deflate allows) decode to
 * {@code invalidSymbol}, which the decoder rejects.
 *
 * One lookup decodes a code of up to {@code primaryBits} bits, which is most
 * symbols; the entry's length is the full code length, so the caller consumes
 * it from the original bits whichever level it came from.
 */
final class Huffman {

	static final int MAX_BITS = 15;

	private Huffman () {}

	/**
	 * Checks code lengths the way zlib does: over-subscribed codes are invalid, and
	 * incomplete ones are only allowed for a single code of length 1 or, if
	 * {@code allowEmpty}, no code at all.
	 */
	static boolean isValid (byte [] lengths, int offset, int n, boolean allowEmpty) {
		int [] counts = new int[MAX_BITS+1];
		for (int i = 0; i < n; i++) counts[lengths[offset+i]]++;
		return isValid(counts, n - counts[0], allowEmpty);
	}

	private static boolean isValid (int [] counts, int codes, boolean allowEmpty) {
		int left = 1;
		for (int len = 1; len <= MAX_BITS; len++) {
			left = (left << 1) - counts[len];
			if (left < 0) return false;
		}
		if (left == 0) return true;
		if (codes == 0) return allowEmpty;
		return codes == 1 && counts[1] == 1;
	}

	/**
	 * Builds the decoding table for the code lengths lengths[offset, offset+n).
	 *
	 * @throws DataFormatError if the lengths aren't a valid code (see isValid)
	 */
	static int [] build (byte [] lengths, int offset, int n, int primaryBits, int invalidSymbol, boolean allowEmpty) throws DataFormatError {
		return build(lengths, offset, n, primaryBits, invalidSymbol, allowEmpty, null);
	}

	/**
	 * Builds the table into reuse if it's large enough, else into a new array. Only
	 * the first size entries are filled and looked up, so a longer reused array's stale
	 * tail is never read.
	 */
	static int [] build (byte [] lengths, int offset, int n, int primaryBits, int invalidSymbol, boolean allowEmpty, int [] reuse) throws DataFormatError {
		int [] counts = new int[MAX_BITS+1];
		for (int i = 0; i < n; i++) counts[lengths[offset+i]]++;
		counts[0] = 0;
		int codes = 0;
		for (int len = 1; len <= MAX_BITS; len++) codes += counts[len];
		if (!isValid(counts, codes, allowEmpty)) throw new DataFormatError("Invalid Huffman code");

		int [] nextCode = new int[MAX_BITS+2];
		int code = 0;
		for (int len = 1; len <= MAX_BITS; len++) {
			code = (code + counts[len-1]) << 1;
			nextCode[len] = code;
		}

		int primarySize = 1 << primaryBits;
		int primaryMask = primarySize - 1;

		// Reversed codes, as deflate packs codes starting from their most significant bit
		int [] reversed = new int[n];
		int [] subBits = null;
		for (int s = 0; s < n; s++) {
			int len = lengths[offset+s];
			if (len == 0) continue;
			reversed[s] = Integer.reverse(nextCode[len]++) >>> (32 - len);
			if (len > primaryBits) {
				if (subBits == null) subBits = new int[primarySize];
				int prefix = reversed[s] & primaryMask;
				subBits[prefix] = Math.max(subBits[prefix], len - primaryBits);
			}
		}

		int size = primarySize;
		int [] subOffset = null;
		if (subBits != null) {
			subOffset = new int[primarySize];
			for (int p = 0; p < primarySize; p++) {
				if (subBits[p] > 0) {
					subOffset[p] = size;
					size += 1 << subBits[p];
				}
			}
		}

		int [] table = reuse != null && reuse.length >= size ? reuse : new int[size];
		java.util.Arrays.fill(table, 0, size, invalidSymbol << 4 | 1);
		for (int s = 0; s < n; s++) {
			int len = lengths[offset+s];
			if (len == 0) continue;
			int entry = s << 4 | len;
			int r = reversed[s];
			if (len <= primaryBits) {
				for (int i = r; i < primarySize; i += 1 << len) table[i] = entry;
			}
			else {
				int prefix = r & primaryMask;
				int base = subOffset[prefix];
				int bits = subBits[prefix];
				table[prefix] = 1 << 31 | base << 4 | bits;
				for (int i = r >>> primaryBits; i < 1 << bits; i += 1 << (len - primaryBits)) table[base + i] = entry;
			}
		}
		return table;
	}

}
