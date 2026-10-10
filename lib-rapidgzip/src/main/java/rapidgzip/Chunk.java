package rapidgzip;

/**
 * The output of decoding the deflate blocks between two block boundaries.
 *
 * Decoded without knowing the preceding 32KB of output, the start of the
 * output is in {@link #markers}: values below 256 are bytes, and 256+i stands
 * for byte i of that unknown window. Once the window is known
 * {@link #resolve(byte[])} turns them into bytes. The rest of the output, once
 * no marker can be referenced any more, is in {@link #bytes}.
 */
final class Chunk {

	/** Bit offset of the first block. */
	final long startBit;
	/** Bit offset of the block boundary where decoding stopped. */
	final long endBit;
	/** True if decoding stopped at the end of the data. */
	final boolean endOfData;

	// Released once resolved, as it's the bulk of the memory in use
	private char [] markers;
	final int markersFrom;
	final int markersTo;

	final byte [] bytes;
	final int bytesFrom;
	final int bytesTo;

	/** Output offset (from the start of this chunk), CRC32 and size of each gzip member ending in this chunk. */
	final long [] memberEnds;
	final int [] memberCrcs;
	final int [] memberSizes;

	Chunk (long startBit, long endBit, boolean endOfData, char [] markers, int markersFrom, int markersTo,
			byte [] bytes, int bytesFrom, int bytesTo, long [] memberEnds, int [] memberCrcs, int [] memberSizes) {
		this.startBit = startBit;
		this.endBit = endBit;
		this.endOfData = endOfData;
		this.markers = markers;
		this.markersFrom = markersFrom;
		this.markersTo = markersTo;
		this.bytes = bytes;
		this.bytesFrom = bytesFrom;
		this.bytesTo = bytesTo;
		this.memberEnds = memberEnds;
		this.memberCrcs = memberCrcs;
		this.memberSizes = memberSizes;
	}

	long length () {
		return (markersTo - markersFrom) + (bytesTo - bytesFrom);
	}

	/** The marker part of the output as bytes, given the up to 32KB window preceding the chunk. */
	byte [] resolve (byte [] window) throws DataFormatError {
		byte [] out = new byte[markersTo - markersFrom];
		resolve(window, markersFrom, markersTo, out, 0);
		return out;
	}

	/** Like resolve, then releases the markers, after which only the resolved output is available. */
	byte [] resolveAndRelease (byte [] window) throws DataFormatError {
		byte [] out = resolve(window);
		markers = null;
		return out;
	}

	/** Resolves markers[from, to) into dest. */
	void resolve (byte [] window, int from, int to, byte [] dest, int offset) throws DataFormatError {
		// Markers index the last 32KB of output before the chunk; a shorter window is right aligned
		int shift = DeflateDecoder.WINDOW - window.length;
		for (int i = from; i < to; i++) {
			char c = markers[i];
			if (c < 256) {
				dest[offset++] = (byte) c;
			}
			else {
				int w = c - 256 - shift;
				if (w < 0) throw new DataFormatError("Invalid distance too far back");
				dest[offset++] = window[w];
			}
		}
	}

	/** The last up to 32KB of output up to the end of this chunk, given the window preceding it. */
	byte [] windowAfter (byte [] window) throws DataFormatError {
		int fromBytes = bytesTo - bytesFrom;
		int fromMarkers = markersTo - markersFrom;
		int size = (int) Math.min(DeflateDecoder.WINDOW, (long) window.length + fromMarkers + fromBytes);
		byte [] next = new byte[size];
		int end = size;
		int n = Math.min(end, fromBytes);
		System.arraycopy(bytes, bytesTo - n, next, end - n, n);
		end -= n;
		n = Math.min(end, fromMarkers);
		resolve(window, markersTo - n, markersTo, next, end - n);
		end -= n;
		System.arraycopy(window, window.length - end, next, 0, end);
		return next;
	}

}
