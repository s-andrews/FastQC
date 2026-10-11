/**
 * Copyright Copyright 2010-17 Simon Andrews
 *
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

import java.io.File;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.charset.Charset;
import java.nio.file.Files;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.compress.compressors.bzip2.BZip2CompressorInputStream;

import rapidgzip.RapidGzipInputStream;
import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Analysis.AnalysisQueue;
import uk.ac.babraham.FastQC.Utilities.MultiMemberGZIPInputStream;
import uk.ac.babraham.FastQC.Utilities.ThreadedInputStream;

public class FastQFile implements SequenceFile {

	private Sequence nextSequence = null;
	private File file;
	private long fileSize = 0;

	private boolean casavaMode = false;
	private boolean nofilter = false;

	// We actually read our final data from this buffered reader
	private ByteLineReader br;
	// Lines are read as bytes; the id line is decoded as the platform charset, as it was read before
	private static final Charset CHARSET = Charset.defaultCharset();
	// We'll keep count of the number of lines read for the error message
	private long lineNumber = 0;

	// We keep the file stream around just so we can see how far through
	// the file we've got.  We don't read from this directly, but it's the
	// only way to access the file pointer.
	private FileInputStream fis;
	// Set when a gzip file is decompressed in parallel. fis is then never read
	// (RapidGzipInputStream opens its own channel), so progress comes from how
	// far decompression got through the compressed file instead.
	private RapidGzipInputStream gzipStream = null;

	private String name;
	private boolean isColorspace = false;

	protected FastQFile(FastQCConfig config,File file) throws SequenceFormatException, IOException {
		this.file = file;
		if (file.getName().startsWith("stdin")) {
			fileSize = Long.MAX_VALUE;
		}
		else {
			fileSize = file.length();
		}
		name = file.getName();

		if (config.casava) {
			casavaMode = true;
			if (config.nofilter) {
				nofilter = true;
			}
		}

		if (!file.getName().startsWith("stdin")) {
			fis = new FileInputStream(file);
		}
				
		if (file.getName().startsWith("stdin")) {
			br = new ByteLineReader(System.in);
		}
		else if (file.getName().toLowerCase().endsWith(".gz") || (Files.probeContentType(file.toPath()) != null && (Files.probeContentType(file.toPath()).equals("application/x-gzip") || Files.probeContentType(file.toPath()).equals("application/gzip")))) {
			// RapidGzipInputStream splits the compressed file into chunks and decodes
			// them speculatively on a thread pool, which needs random access to the
			// file and spare cores. Decompression was the reader's bottleneck: wall
			// -5% on AVITI 2M reads, and it lets more processor threads pay off.
			// The sequential paths use a 64 KB input buffer: the 512 byte default means
			// one native read per 512 compressed bytes (wall -2% on 2M AVITI reads).
			if (!Files.isRegularFile(file.toPath())) {
				// Eg a named pipe (process substitution), which can only be read sequentially
				br = new ByteLineReader(decompressor(config, new MultiMemberGZIPInputStream(fis, 64 * 1024)));
			}
			else if (config.threads != null && config.threads <= 1) {
				// Single-threaded as requested: decompress inline on the reading thread
				br = new ByteLineReader(new MultiMemberGZIPInputStream(fis, 64 * 1024));
			}
			else {
				// Decompress on several threads; the decoded output is still read in
				// order, so the sequences and the report are identical. It reads the
				// file itself, so fis isn't needed (and can't leak if this throws).
				fis.close();
				fis = null;
				gzipStream = new RapidGzipInputStream(file.toPath(), gzipThreads(config));
				br = new ByteLineReader(gzipStream);
			}
		} 
		else if (file.getName().toLowerCase().endsWith(".bz2")) {
			br = new ByteLineReader(decompressor(config, new BZip2CompressorInputStream(fis,true)));
		} 

		else {
			br = new ByteLineReader(fis);
		}
		readNext();
	}

	/**
	 * Threads to decompress a gzip file with: half of this file's share of the
	 * thread budget if one was given (-t is shared by the files processed in
	 * parallel), else 4, as beyond that decompression is no longer the bottleneck
	 * (decompressing AVITI alone: 3.6 s with GZIPInputStream, 1.9 s with 4 threads,
	 * 1.4 s with 8 for more CPU). These threads come on top of the per-file
	 * reader and processor threads sized by AnalysisQueue.
	 */
	private static int gzipThreads(FastQCConfig config) {
		if (config.threads != null) return Math.max(1, config.threads / AnalysisQueue.getInstance().getFilesInParallel() / 2);
		return 4;
	}

	/**
	 * Decompression is as expensive as parsing, so unless we're limited to a
	 * single thread run it on its own thread, overlapping with the parsing
	 * (wall -10% on 2M AVITI reads). With -t 1 the stream is read directly
	 * on the calling thread. The extra thread isn't counted in -t.
	 * Note that the progress, taken from the position of fis, runs ahead of
	 * the parsing by the data buffered in ThreadedInputStream.
	 */
	private InputStream decompressor(FastQCConfig config, InputStream in) {
		if (config.threads != null && config.threads <= 1) return in;
		return new ThreadedInputStream(in, "fastqc-decompress");
	}

	public String name() {
		return name;
	}

	public int getPercentComplete() {
		if (! hasNext()) return 100;
		if (file.getName().startsWith("stdin")) {
			return 0;
		}
		if (gzipStream != null) {
			// Backed by a volatile field the consumer advances, so a cheap, safe poll
			return (int) (((double)gzipStream.compressedPosition()/ fileSize)*100);
		}
		try {
			int percent = (int) (((double)fis.getChannel().position()/ fileSize)*100);
			return percent;
		} 
		catch (IOException e) {
			e.printStackTrace();
		}
		return 0;
	}

	public boolean isColorspace() {
		return isColorspace;
	}

	public void setIsColorspace(boolean isColorspace) {
		this.isColorspace = isColorspace;
	}

	public boolean hasNext() {
		return nextSequence != null;
	}

	public Sequence next() throws SequenceFormatException {
		Sequence seq = nextSequence;
		readNext();
		return seq;
	}

	private void readNext() throws SequenceFormatException {
		try {
			readNextRecord();
		}
		catch (SequenceFormatException | RuntimeException e) {
			// The file can't be read any further and nothing else closes it: release
			// the file and the decompression thread, which would otherwise stay
			// blocked for the rest of the run with its buffered chunks
			try {
				br.close();
				if (fis != null) fis.close();
			}
			catch (IOException ignored) {}
			throw e;
		}
	}

	private void readNextRecord() throws SequenceFormatException {
		try {
			// First line should be the id

			// We might have blank lines between entries or at the end
			// so allow for this
			byte [] idLine;

			while (true) {
				idLine = br.readLine();
				++lineNumber;

				if (idLine == null) {
					nextSequence = null;
					br.close();
					if (fis != null) {
						fis.close();
					}
					return;
				}
				if (idLine.length==0) {
					continue;
				}				

				break;
			}

			String id = new String(idLine, CHARSET);

			if (!id.startsWith("@")) {
				nextSequence = null;
				throw new SequenceFormatException("ID line didn't start with '@' at line "+lineNumber);
			}

			byte [] seq;
			byte [] midLine;
			byte [] quality;
			try {
				// Then the sequence
				seq = br.readLine();
				if (seq == null) throw new IOException("No more data, expected sequence, at line "+lineNumber);
				++lineNumber;
				// Then another id which we don't need
				midLine = br.readLine();
				if (midLine == null) throw new IOException("No more data, expected midline, at line "+lineNumber);
				++lineNumber;
				if (midLine.length == 0 || midLine[0] != '+') {
					throw new SequenceFormatException("Midline '"+new String(midLine, CHARSET)+"' didn't start with '+' at "+lineNumber);
				}
				// Then the quality string
				quality = br.readLine();
				if (quality == null) throw new IOException("No more data, expected quality, at line "+lineNumber);
				++lineNumber;
			}
			catch (IOException ioe) {
				throw new SequenceFormatException("Ran out of data in the middle of a fastq entry.  Your file is probably truncated");
			}


			// We only check for colourspace on the first entry.  After that we assume
			// the rest of the file is the same.  For the first entry the nextSequence
			// will be null, but we'll have real data in seq
			if (nextSequence == null && seq != null) {
				checkColorspace(new String(seq, CHARSET));
			}

			if (isColorspace()) {
				String seqString = new String(seq, CHARSET);
				nextSequence = new Sequence(this,convertColorspaceToBases(seqString.toUpperCase()), seqString.toUpperCase(), new String(quality, CHARSET), id);
			} 
			else {
				// The lines are fresh copies, so they can be upper cased in place and handed over
				nextSequence = new Sequence(this, toUpperCase(seq), quality, id);
			}

			// If we're running in --casava mode then we will flag any sequences which
			// are marked as being filtered.
			if (casavaMode  && !nofilter) {

				// This is the test illumina suggest, but it's a bit flakey, and I'm not
				// sure it's not going to catch things it shouldn't.
				if (id.indexOf(":Y:") > 0) {
					nextSequence.setIsFiltered(true);
				}
			}




		} 
		// We could have other errors like zip format exceptions which we need to catch
		catch (IOException ioe) {
			throw new SequenceFormatException(ioe.getLocalizedMessage());
		}
	}

	/**
	 * Upper cases ASCII letters in place. Same as String.toUpperCase for the
	 * ASCII bases of a valid FASTQ file (and, unlike it, independent of locale).
	 */
	private static byte [] toUpperCase(byte [] seq) {
		for (int i=0;i<seq.length;i++) {
			byte b = seq[i];
			if (b >= 'a' && b <= 'z') seq[i] = (byte) (b - 32);
		}
		return seq;
	}

	private void checkColorspace(String seq) {
		// Some basecalled files can be all dots, which leads to them
		// being identified as colorspace data. This check should find
		// only true colorspace files.
		String regex = "^[GATCNgatcn][\\.0123456]+$";
		Pattern pattern = Pattern.compile(regex);
		Matcher matcher = pattern.matcher(seq);
		if (matcher.find()) {
			isColorspace = true;
		} else {
			isColorspace = false;
		}
	}

	private String convertColorspaceToBases(String s) {

		char[] cs = s.toUpperCase().toCharArray();

		// We've had a crash report where a file contained a zero length
		// colorspace entry.  This is completely invalid, but we should
		// handle it anyway.
		if (cs.length == 0) {
			return "";
		}

		char[] bp = new char[cs.length - 1];

		char refBase;

		for (int i = 1; i < cs.length; i++) {
			if (i == 1) {
				refBase = cs[i - 1];
			} else {
				refBase = bp[i - 2];
			}
			if (!(refBase == 'G' || refBase == 'A' || refBase == 'T' || refBase == 'C')) {
				throw new IllegalArgumentException("Colourspace sequence data should always start with a real DNA letter.  Line '"+s+"' started with " + refBase
						+ " at position " + i);
			}
			switch (cs[i]) {
			case ('0'):
				switch (refBase) {
				case ('G'):
					bp[i - 1] = 'G';
				break;
				case ('A'):
					bp[i - 1] = 'A';
				break;
				case ('T'):
					bp[i - 1] = 'T';
				break;
				case ('C'):
					bp[i - 1] = 'C';
				break;
				}
			break;
			case ('1'):
				switch (refBase) {
				case ('G'):
					bp[i - 1] = 'T';
				break;
				case ('A'):
					bp[i - 1] = 'C';
				break;
				case ('T'):
					bp[i - 1] = 'G';
				break;
				case ('C'):
					bp[i - 1] = 'A';
				break;
				}
			break;

			case ('2'):
				switch (refBase) {
				case ('G'):
					bp[i - 1] = 'A';
				break;
				case ('A'):
					bp[i - 1] = 'G';
				break;
				case ('T'):
					bp[i - 1] = 'C';
				break;
				case ('C'):
					bp[i - 1] = 'T';
				break;
				}
			break;

			case ('3'):
				switch (refBase) {
				case ('G'):
					bp[i - 1] = 'C';
				break;
				case ('A'):
					bp[i - 1] = 'T';
				break;
				case ('T'):
					bp[i - 1] = 'A';
				break;
				case ('C'):
					bp[i - 1] = 'G';
				break;
				}
			break;

			case ('.'):
			case ('4'):
			case ('5'):
			case ('6'):
				for (; i < cs.length; i++) {
					bp[i - 1] = 'N';
				}
			break;
			default:
				throw new IllegalArgumentException("Unexpected cs char "
						+ cs[i]);
			}
		}

		return new String(bp);
	}

	public void remove() {
		// No action here
	}

	public File getFile() {
		return file;
	}

}
