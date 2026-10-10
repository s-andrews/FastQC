/**
 *    This file is part of FastQC.
 *
 *    FastQC is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    FastQC is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with FastQC; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
package uk.ac.babraham.FastQC.Sequence;

import java.io.IOException;
import java.io.InputStream;
import java.lang.invoke.MethodHandles;
import java.lang.invoke.VarHandle;
import java.nio.ByteOrder;
import java.util.Arrays;

/**
 * Reads lines as bytes, without decoding them into characters. Like
 * BufferedReader.readLine, a line ends at '\n', '\r' or "\r\n".
 */
class ByteLineReader {

	private static final VarHandle LONG_LE = MethodHandles.byteArrayViewVarHandle(long[].class, ByteOrder.LITTLE_ENDIAN);
	private static final long ONES = 0x0101010101010101L;
	private static final long HIGHS = 0x8080808080808080L;
	private static final long LFS = ONES * '\n';
	private static final long CRS = ONES * '\r';

	private final InputStream in;
	private byte [] buffer = new byte[1 << 16];
	private int pos = 0;
	private int limit = 0;
	private boolean endOfStream = false;
	// The last line ended with '\r', so a following '\n' is part of that line end
	private boolean skipLF = false;

	ByteLineReader (InputStream in) {
		this.in = in;
	}

	/** The next line without its terminator, or null at the end of the stream. */
	byte [] readLine () throws IOException {
		int scanned = 0;
		while (true) {
			if (skipLF) {
				if (pos == limit && !fill()) return null;
				if (buffer[pos] == '\n') pos++;
				skipLF = false;
			}
			int end = findLineEnd(pos + scanned);
			if (end >= 0) {
				byte [] line = Arrays.copyOfRange(buffer, pos, end);
				if (buffer[end] == '\r') skipLF = true;
				pos = end + 1;
				return line;
			}
			scanned = limit - pos;
			if (!fill()) {
				if (pos == limit) return null;
				// Last line without a terminator
				byte [] line = Arrays.copyOfRange(buffer, pos, limit);
				pos = limit;
				return line;
			}
		}
	}

	/** Index of the first '\n' or '\r' from start, or -1. */
	private int findLineEnd (int start) {
		int i = start;
		// 8 bytes at a time: a byte equal to LF or CR becomes zero after the XOR
		for (; i + 8 <= limit; i += 8) {
			long word = (long) LONG_LE.get(buffer, i);
			long lf = word ^ LFS;
			long cr = word ^ CRS;
			long found = ((lf - ONES) & ~lf | (cr - ONES) & ~cr) & HIGHS;
			if (found != 0) return i + (Long.numberOfTrailingZeros(found) >>> 3);
		}
		for (; i < limit; i++) {
			if (buffer[i] == '\n' || buffer[i] == '\r') return i;
		}
		return -1;
	}

	/** Reads more data, keeping the unread part. Returns false at the end of the stream. */
	private boolean fill () throws IOException {
		if (endOfStream) return false;
		if (pos > 0) {
			System.arraycopy(buffer, pos, buffer, 0, limit - pos);
			limit -= pos;
			pos = 0;
		}
		if (limit == buffer.length) buffer = Arrays.copyOf(buffer, buffer.length * 2);
		int n = in.read(buffer, limit, buffer.length - limit);
		if (n < 0) {
			endOfStream = true;
			return false;
		}
		limit += n;
		return true;
	}

	void close () throws IOException {
		in.close();
	}

}
