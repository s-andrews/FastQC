package rapidgzip;

import java.util.zip.ZipException;

/** Invalid deflate or gzip data. */
final class DataFormatError extends ZipException {

	private static final long serialVersionUID = 1L;

	DataFormatError (String message) {
		super(message);
	}

}
