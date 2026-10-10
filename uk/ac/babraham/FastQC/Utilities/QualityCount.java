/**
 * Copyright Copyright 2013-17 Simon Andrews
 *
 *    This file is part of FastQC.
 *
 *    Conclave is free software; you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation; either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    SeqMonk is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with Conclave; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

package uk.ac.babraham.FastQC.Utilities;

public class QualityCount {

	/*
	 * So I'm on my third go at writing this.  I've now tried an all
	 * primitive version of this class so that we don't have to do 
	 * hash lookups which require a conversion from chr to Character.
	 * We should also be safe with 150 slots which will give us up to
	 * Phred 86 with a 64 offset, which should be plenty.
	 */
	
	/** Quality characters must be below this. */
	public static final int MAX_CHAR = 150;

	// Counts of the characters from first on, covering only the range seen so far:
	// there are thousands of these (per tile and position) and real qualities only
	// span a few dozen characters
	private long [] actualCounts = new long[0];
	private int first = 0;
	
	private long totalCounts = 0;

	public void addValue(char c) {
		totalCounts++;
		
		if ((int)c >= MAX_CHAR) {
			throw new IllegalArgumentException("Got character "+c+" as a quality value which has ASCII"+((int)c)+" which is higher than we can cope with");
		}
		include(c);
		actualCounts[c - first]++;
	}
	
	/**
	 * Adds count occurrences of c, which must be below MAX_CHAR (not checked:
	 * callers count into flat tables first and check characters there).
	 */
	public void addValue(char c, long count) {
		totalCounts += count;
		include(c);
		actualCounts[c - first] += count;
	}

	/** Extends the counted range to include c. */
	private void include (int c) {
		if (actualCounts.length == 0) {
			first = c;
			actualCounts = new long[1];
		}
		else if (c < first) {
			long [] counts = new long[actualCounts.length + first - c];
			System.arraycopy(actualCounts, 0, counts, first - c, actualCounts.length);
			actualCounts = counts;
			first = c;
		}
		else if (c >= first + actualCounts.length) {
			actualCounts = java.util.Arrays.copyOf(actualCounts, c - first + 1);
		}
	}

	private long count (int c) {
		int i = c - first;
		return i >= 0 && i < actualCounts.length ? actualCounts[i] : 0;
	}
	
	public long getTotalCount () {
		return totalCounts;
	}
	
	public char getMinChar () {
		
		for (int i=0;i<MAX_CHAR;i++) {
			if (count(i)>0) return (char)i;
		}
		
		return (char)1000;
	}
	
	public char getMaxChar () {
		for (int i=MAX_CHAR-1;i>=0;i--) {
			if (count(i)>0) return (char)i;
		}
		
		return (char)1000;

	}
			
	public double getMean (int offset) {
		long total = 0;
		long count = 0;
	
		for (int i=offset;i<MAX_CHAR;i++) {
			total += count(i)*(i-offset);
			count += count(i);
		}
		
		return ((double)total)/count;
	}
	
	public double getPercentile (int offset, int percentile) {

		long total = totalCounts;
		
		total *= percentile;
		total /= 100;
		
		long count = 0;
		for (int i=offset;i<MAX_CHAR;i++) {
			count += count(i);
			if (count >=total) {
				return((char)(i-offset));
			}
		}
		
		return -1;
		
	}
	
}
