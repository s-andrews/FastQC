package rapidgzip;

import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Reuses the large output buffers of chunks: each chunk otherwise allocates
 * several MB, which keeps the garbage collector busy. Buffers handed out may
 * hold stale data and be larger than asked, so users must keep track of the
 * range they wrote.
 *
 * Shared by the decoding and resolving threads of the pool and the reading
 * thread, hence the concurrent queues. A buffer must only be released by its
 * last user once nothing reads it any more: arrays a DeflateDecoder outgrew,
 * a chunk's markers once resolved, and output segments once read.
 * At most max buffers of each type are kept (one per chunk in flight).
 */
final class Buffers {

	// Smaller buffers aren't worth keeping
	private static final int MIN_SIZE = 1 << 16;

	private final int max;
	private final ConcurrentLinkedQueue<char []> chars = new ConcurrentLinkedQueue<char []>();
	private final ConcurrentLinkedQueue<byte []> bytes = new ConcurrentLinkedQueue<byte []>();

	Buffers (int max) {
		this.max = max;
	}

	/** A buffer of at least size values; pooled ones too small are dropped. */
	char [] chars (int size) {
		char [] c;
		while ((c = chars.poll()) != null) {
			if (c.length >= size) return c;
		}
		return new char[size];
	}

	byte [] bytes (int size) {
		byte [] b;
		while ((b = bytes.poll()) != null) {
			if (b.length >= size) return b;
		}
		return new byte[size];
	}

	// The size check and the add aren't atomic, so the pool may go slightly over max
	void release (char [] c) {
		if (c != null && c.length >= MIN_SIZE && chars.size() < max) chars.add(c);
	}

	void release (byte [] b) {
		if (b != null && b.length >= MIN_SIZE && bytes.size() < max) bytes.add(b);
	}

}
