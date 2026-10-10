package rapidgzip;

import java.util.concurrent.ConcurrentLinkedQueue;

/**
 * Reuses the large output buffers of chunks: each chunk otherwise allocates
 * several MB, which keeps the garbage collector busy. Buffers handed out may
 * hold stale data.
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

	void release (char [] c) {
		if (c != null && c.length >= MIN_SIZE && chars.size() < max) chars.add(c);
	}

	void release (byte [] b) {
		if (b != null && b.length >= MIN_SIZE && bytes.size() < max) bytes.add(b);
	}

}
