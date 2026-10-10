package rapidgzip;

/**
 * The output of decoding the deflate blocks between two block boundaries.
 *
 * Decoded without knowing the preceding 32KB of output, the start of the
 * output is in {@link #markers}: values below 256 are bytes, and 256+i stands
 * for byte i of that unknown window. Once the window is known
 * {@link #resolve(byte[])} turns them into bytes. The rest of the output, once
 * no marker can be referenced any more, is in {@link #bytes}.
 *
 * Built by a worker, then read by the reader thread (windowAfter) and resolved
 * on a worker; fields are final apart from markers, which are released after
 * resolving as they are the bulk of the memory in flight (2 bytes per output
 * byte). Resolving is a plain lookup, so it gives exactly the bytes a decoder
 * knowing the window would have produced.
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

	int markersLength () {
		return markersTo - markersFrom;
	}

	/**
	 * The marker part of the output as bytes, given the up to 32KB window
	 * preceding the chunk, in the first markersLength() bytes of a buffer. The
	 * markers are released, after which only the resolved output is available.
	 */
	byte [] resolveAndRelease (byte [] window, Buffers buffers) throws DataFormatError {
		byte [] out = buffers.bytes(markersLength());
		resolve(window, markersFrom, markersTo, out, 0);
		buffers.release(markers);
		markers = null;
		return out;
	}

	/**
	 * Resolves markers[from, to) into dest. Called twice per speculative chunk (for the
	 * next window, then for the whole marker output on the pool); both only read markers.
	 */
	void resolve (byte [] window, int from, int to, byte [] dest, int offset) throws DataFormatError {
		// A table from any 16 bit value to its byte: values are as often bytes as markers, so
		// a lookup avoids a mispredicted branch for each. Markers index the last 32KB of output
		// before the chunk; a shorter window is right aligned. A marker before the start of a
		// short window is a distance a sequential decoder rejects too.
		// The 33KB table is built per call, which is small next to a chunk's output.
		byte [] values = new byte[256 + DeflateDecoder.WINDOW];
		for (int i = 0; i < 256; i++) values[i] = (byte) i;
		int shift = DeflateDecoder.WINDOW - window.length;
		System.arraycopy(window, 0, values, 256 + shift, window.length);
		// Markers are below 256 + WINDOW by construction
		char [] in = markers;
		for (int i = from; i < to; i++) {
			dest[offset++] = values[in[i]];
		}
		if (shift > 0) {
			// Markers before the start of the data are invalid distances. Checked in a separate
			// pass to keep the loop above branch free; it only runs while less than 32KB of
			// output precedes the chunk. The bytes written above for them are then never used.
			for (int i = from; i < to; i++) {
				char c = markers[i];
				if (c >= 256 && c < 256 + shift) throw new DataFormatError("Invalid distance too far back");
			}
		}
	}

	/**
	 * The last up to 32KB of output up to the end of this chunk, given the window preceding it.
	 * Resolves only the markers it needs, so the next chunk's window is known without waiting
	 * for the whole chunk to be resolved: this is what lets resolving run in parallel.
	 */
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
