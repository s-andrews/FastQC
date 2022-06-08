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
package uk.ac.babraham.FastQC.Sequence.QualityEncoding;

public class PhredEncoding {

	private String name;
	private int offset;
	
	private static final int SANGER_ENCODING_OFFSET = 33;
	private static final int ILLUMINA_1_3_ENCODING_OFFSET = 64;
	
	public static PhredEncoding getFastQEncodingOffset (char lowestChar) {
		if (lowestChar < 33) {
			throw new IllegalArgumentException("No known encodings with chars < 33 (Yours was '"+lowestChar+"' with value "+(int)lowestChar+")");
		}
		else if (lowestChar < 64) {
			return new PhredEncoding("Sanger / Illumina 1.9", SANGER_ENCODING_OFFSET);
		}
		
		// There are potentially two encodings using an offset of 64.  Illumina
		// v1.3 allowed quality values of 1, whereas from v1.5 onwards the lowest
		// value allowed was 2.  If we guess wrong between these two then it's not
		// the end of the world since they use the same offset.
		else if (lowestChar == ILLUMINA_1_3_ENCODING_OFFSET+1) {
			return new PhredEncoding("Illumina 1.3", ILLUMINA_1_3_ENCODING_OFFSET);			
		}
		else if (lowestChar <= 126) {
			return new PhredEncoding("Illumina 1.5", ILLUMINA_1_3_ENCODING_OFFSET);
		}
		throw new IllegalArgumentException("No known encodings with chars > 126 (Yours was "+lowestChar+" with value "+(int)lowestChar+")");
	}
	
	public static double convertSangerPhredToProbability (int phred) {
		return Math.pow(10,phred/-10d);
	}
	
	public static double convertOldIlluminaPhredToProbability (int phred) {
		return Math.pow(10, ((double)phred/(phred+1))/-10d);
	}
	
	public static int convertProbabilityToSangerPhred (double p) {
		return (int)Math.round(-10d*Math.log10(p));
	}

	public static int convertProbabilityToOldIlluminaPhred (double p) {
		return (int)Math.round(-10d*Math.log10(p/1-p));
	}

	private PhredEncoding (String name, int offset) {
		this.name = name;
		this.offset = offset;
	}
	
	public String name () {
		return name;
	}
	
	public String toString () {
		return name();
	}
	
	public int offset () {
		return offset;
	}
	
	
	public static void main (String [] args) {
		double p = 0.4;
		
		System.out.println("Sanger phred for p="+p+" is "+convertProbabilityToSangerPhred(p));
		
		int phred=4;
		System.out.println("P value for Sanger phred="+phred+" is "+convertSangerPhredToProbability(phred));

	}
}
