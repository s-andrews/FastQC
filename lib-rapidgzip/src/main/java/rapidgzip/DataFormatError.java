package rapidgzip;

import java.util.zip.ZipException;

/**
 * Invalid deflate or gzip data. Thrown when decoding from a guessed block start
 * means the guess was wrong; a ZipException, like GZIPInputStream's errors.
 */
final class DataFormatError extends ZipException {

	private static final long serialVersionUID = 1L;

	DataFormatError (String message) {
		super(message);
	}

}
