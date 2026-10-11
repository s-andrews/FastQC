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
 * <p>
 * Inflating took about a quarter of the CPU time on the reader thread, which
 * was the critical path: overlapping it with parsing cut wall time by ~10%
 * (2M AVITI reads, 3.61 s to 3.25 s). The bytes delivered are exactly those
 * of the wrapped stream, in order.
 * <p>
 * Threading: one producer thread (created here) reads the wrapped stream and
 * closes it; the read methods and close() must be called from a single
 * consumer thread. At most QUEUED_CHUNKS chunks plus one being read and one
 * being consumed are held in memory (~1.5 MB).
 */
public class ThreadedInputStream extends InputStream {

	// Large chunks keep the hand-off cost (one queue operation per chunk)
	// negligible; a few queued chunks absorb jitter between the two threads.
	private static final int CHUNK_SIZE = 256 * 1024;
	private static final int QUEUED_CHUNKS = 4;

	// End-of-stream marker, compared by identity. It is queued after the last
	// chunk both on a normal end and when the producer hit an IOException.
	private static final byte [] EOF = new byte[0];

	private final BlockingQueue<byte []> queue = new ArrayBlockingQueue<byte []>(QUEUED_CHUNKS);
	private final Thread producer;
	// Set by the producer before it queues EOF, so the consumer sees it
	// once it has taken EOF (the queue hand-off also orders the two).
	private volatile IOException error = null;

	// Consumer-side state, only touched by the consumer thread.
	private byte [] chunk = null;
	private int chunkLength = 0;
	private int position = 0;
	private boolean finished = false;

	public ThreadedInputStream (final InputStream in, String name) {
		producer = new Thread(new Runnable() {
			public void run() {
				try {
					while (true) {
						// readNBytes blocks until the chunk is full or the stream ends,
						// so only the last chunk can be short; empty means end of stream.
						byte [] buffer = in.readNBytes(CHUNK_SIZE);
						if (buffer.length == 0) break;
						queue.put(buffer);
					}
				}
				catch (IOException e) {
					error = e;
				}
				catch (InterruptedException e) {
					// The reader was closed: nobody will take EOF, so don't queue it
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
		// Daemon, so a stream that is never closed can't keep the JVM alive.
		producer.setDaemon(true);
		producer.start();
	}

	/**
	 * Takes the next chunk from the producer, blocking until one is ready.
	 * Returns false at the end of the stream; a producer IOException is
	 * thrown here once, after all the data read before it.
	 */
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

	// Returns at most what's left of the current chunk, which is allowed by
	// the InputStream contract; callers like InputStreamReader loop as needed.
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

	/**
	 * Stops the producer: interrupting it ends a blocked put(), and the
	 * producer closes the wrapped stream itself, so close() never races with
	 * a read in progress on it.
	 */
	public void close () {
		finished = true;
		producer.interrupt();
	}

}
