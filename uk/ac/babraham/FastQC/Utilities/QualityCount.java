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
	
	private long [] actualCounts = new long[150];
	
	private long totalCounts = 0;

	public void addValue(char c) {
		totalCounts++;
		
		if ((int)c >= actualCounts.length) {
			throw new IllegalArgumentException("Got character "+c+" as a quality value which has ASCII"+((int)c)+" which is higher than we can cope with");
		}
		actualCounts[(int)c]++;
	}
	
	public long getTotalCount () {
		return totalCounts;
	}
	
	public char getMinChar () {
		
		for (int i=0;i<actualCounts.length;i++) {
			if (actualCounts[i]>0) return (char)i;
		}
		
		return (char)1000;
	}
	
	public char getMaxChar () {
		for (int i=actualCounts.length-1;i>=0;i--) {
			if (actualCounts[i]>0) return (char)i;
		}
		
		return (char)1000;

	}
			
	public double getMean (int offset) {
		long total = 0;
		long count = 0;
	
		for (int i=offset;i<actualCounts.length;i++) {
			total += actualCounts[i]*(i-offset);
			count += actualCounts[i];
		}
		
		return ((double)total)/count;
	}
	
	public double getPercentile (int offset, int percentile) {

		long total = totalCounts;
		
		total *= percentile;
		total /= 100;
		
		long count = 0;
		for (int i=offset;i<actualCounts.length;i++) {
			count += actualCounts[i];
			if (count >=total) {
				return((char)(i-offset));
			}
		}
		
		return -1;
		
	}
	
}
