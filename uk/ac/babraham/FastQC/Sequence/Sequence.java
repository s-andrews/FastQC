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

import java.nio.charset.StandardCharsets;

/**
 * One read. Bases and qualities are held as byte arrays, which the per-base
 * modules read directly: indexing a byte[] is cheaper than String.charAt and
 * FASTQ records need no decoding into Strings at all (see ByteLineReader).
 * Strings are only built when a module asks for them. A Sequence is not
 * modified once created, so modules on several threads can share it.
 */
public class Sequence {

	// Bases (upper case) and qualities, one byte per position
	private final byte [] bases;
	private final byte [] qualities;
	// The strings this was created from, if it was
	private final String sequence;
	private final String quality;
	private final String id;
	private SequenceFile file;
	private String colorspace;
	private boolean isFiltered;
	
	public Sequence (SequenceFile file,String sequence, String quality, String id) {
		this.id = id;
		this.file = file;
		this.sequence = sequence.toUpperCase();
		this.quality = quality;
		this.bases = toBytes(this.sequence);
		this.qualities = toBytes(quality);
		this.colorspace = null;
		this.isFiltered = false;
	}
	
	public Sequence (SequenceFile file,String sequence, String colorspace, String quality, String id) {
		this.id = id;
		this.file = file;
		this.sequence = sequence;
		this.quality = quality;
		this.bases = toBytes(sequence);
		this.qualities = toBytes(quality);
		this.colorspace = colorspace;
	}

	/**
	 * Creates a sequence from bytes, one per position, which it then owns.
	 * Bases must already be upper case.
	 */
	public Sequence (SequenceFile file, byte [] bases, byte [] qualities, String id) {
		this.id = id;
		this.file = file;
		this.bases = bases;
		this.qualities = qualities;
		this.sequence = null;
		this.quality = null;
		this.colorspace = null;
		this.isFiltered = false;
	}

	// Characters beyond Latin-1 become '?', they're only kept in the strings
	private static byte [] toBytes (String s) {
		return s.getBytes(StandardCharsets.ISO_8859_1);
	}
	
	public void setIsFiltered (boolean isFiltered) {
		this.isFiltered = isFiltered;
	}
	
	public boolean isFiltered () {
		return isFiltered;
	}
	
	public SequenceFile file () {
		return file;
	}

	public int length () {
		return bases.length;
	}

	/** The bases, upper case, one byte each. Must not be modified. */
	public byte [] getSequenceBytes () {
		return bases;
	}

	/** The quality characters, one byte each. Must not be modified. */
	public byte [] getQualityBytes () {
		return qualities;
	}
	
	/** Built on every call when this was created from bytes; prefer getSequenceBytes in per-read loops. */
	public String getSequence () {
		return sequence != null ? sequence : new String(bases, StandardCharsets.ISO_8859_1);
	}
	
	public String getColorspace () {
		return colorspace;
	}
	
	public String getQualityString () {
		return quality != null ? quality : new String(qualities, StandardCharsets.ISO_8859_1);
	}
	
	public String getID () {
		return id;
	}
	
}
