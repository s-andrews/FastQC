package rapidgzip;

import java.io.EOFException;
import java.io.IOException;
import java.util.zip.CRC32;

/** Gzip member header (RFC 1952). */
final class GzipHeader {

	private static final int FHCRC = 2;
	private static final int FEXTRA = 4;
	private static final int FNAME = 8;
	private static final int FCOMMENT = 16;

	private GzipHeader () {}

	/**
	 * Skips a gzip header at the current, byte aligned, position. Returns false if
	 * there's no valid header there, including a header CRC mismatch.
	 * Truncation also returns false: the reader's zero padding ends FNAME/FCOMMENT,
	 * then pastEnd() is checked. After the first member, false means the end of
	 * the data, as GZIPInputStream ignores what follows.
	 */
	static boolean skip (BitReader in) throws IOException {
		try {
			CRC32 crc = new CRC32();
			if (read(in, crc) != 0x1f || read(in, crc) != 0x8b) return false;
			if (read(in, crc) != 8) return false;
			int flags = read(in, crc);
			// mtime, extra flags and OS
			for (int i = 0; i < 6; i++) read(in, crc);
			if ((flags & FEXTRA) != 0) {
				int length = read(in, crc) | read(in, crc) << 8;
				for (int i = 0; i < length; i++) read(in, crc);
			}
			if ((flags & FNAME) != 0) {
				while (read(in, crc) != 0) {}
			}
			if ((flags & FCOMMENT) != 0) {
				while (read(in, crc) != 0) {}
			}
			if ((flags & FHCRC) != 0) {
				int expected = (int) crc.getValue() & 0xFFFF;
				if ((in.readByte() | in.readByte() << 8) != expected) return false;
			}
			return !in.pastEnd();
		}
		catch (EOFException e) {
			return false;
		}
	}

	private static int read (BitReader in, CRC32 crc) throws IOException {
		int b = in.readByte();
		crc.update(b);
		return b;
	}

}
