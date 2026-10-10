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
package uk.ac.babraham.FastQC.Utilities;

import java.io.IOException;
import java.io.InputStream;
import java.io.InterruptedIOException;
import java.util.concurrent.ArrayBlockingQueue;
import java.util.concurrent.BlockingQueue;

/**
 * Reads an input stream ahead on a background thread, so that expensive
 * decoding (eg gzip or bzip2 decompression) overlaps with whatever consumes
 * the data. Exceptions from the underlying stream are rethrown to the reader.
 */
public class ThreadedInputStream extends InputStream {

	private static final int CHUNK_SIZE = 256 * 1024;
	private static final int QUEUED_CHUNKS = 4;

	private static final byte [] EOF = new byte[0];

	private final BlockingQueue<byte []> queue = new ArrayBlockingQueue<byte []>(QUEUED_CHUNKS);
	private final Thread producer;
	private volatile IOException error = null;

	private byte [] chunk = null;
	private int chunkLength = 0;
	private int position = 0;
	private boolean finished = false;

	public ThreadedInputStream (final InputStream in, String name) {
		producer = new Thread(new Runnable() {
			public void run() {
				try {
					while (true) {
						byte [] buffer = in.readNBytes(CHUNK_SIZE);
						if (buffer.length == 0) break;
						queue.put(buffer);
					}
				}
				catch (IOException e) {
					error = e;
				}
				catch (InterruptedException e) {
					// The reader was closed
					return;
				}
				finally {
					try {
						in.close();
					}
					catch (IOException e) {}
				}
				try {
					queue.put(EOF);
				}
				catch (InterruptedException e) {}
			}
		}, name);
		producer.setDaemon(true);
		producer.start();
	}

	private boolean nextChunk () throws IOException {
		if (finished) return false;
		try {
			chunk = queue.take();
		}
		catch (InterruptedException e) {
			throw new InterruptedIOException();
		}
		if (chunk == EOF) {
			finished = true;
			if (error != null) throw error;
			return false;
		}
		chunkLength = chunk.length;
		position = 0;
		return true;
	}

	public int read () throws IOException {
		if (position == chunkLength && !nextChunk()) return -1;
		return chunk[position++] & 0xFF;
	}

	public int read (byte [] b, int off, int len) throws IOException {
		if (len == 0) return 0;
		if (position == chunkLength && !nextChunk()) return -1;
		int n = Math.min(len, chunkLength - position);
		System.arraycopy(chunk, position, b, off, n);
		position += n;
		return n;
	}

	public int available () {
		return chunkLength - position;
	}

	public void close () {
		finished = true;
		producer.interrupt();
	}

}
