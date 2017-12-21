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
package uk.ac.babraham.FastQC.Statistics;


/**
 * A Class to calculate the Pearson Correlation.
 */
public class PearsonCorrelation {

	
	/**
	 * Calculate correlation.
	 * 
	 * @param data1 the first dataset
	 * @param data2 the second dataset
	 * @return the Pearson r-value
	 * @throws SeqMonkException if the two datasets don't have the same number of points in them.
	 */
	public static float calculateCorrelation (long [] data1, long [] data2) {

		float [] d1 = new float[data1.length];
		float [] d2 = new float[data2.length];
		for (int i=0;i<data1.length;i++)d1[i] = data1[i];
		for (int i=0;i<data2.length;i++)d2[i] = data2[i];
		
		return calculateCorrelation(d1, d2);
		
	}
	
	/**
	 * Calculate correlation.
	 * 
	 * @param data1 the first dataset
	 * @param data2 the second dataset
	 * @return the Pearson r-value
	 * @throws SeqMonkException if the two datasets don't have the same number of points in them.
	 */
	public static float calculateCorrelation (long [] data1, long [] data2, int offset) {

		float [] d1 = new float[data1.length-offset];
		float [] d2 = new float[data2.length-offset];
		for (int i=0;i<d1.length;i++)d1[i] = data1[i];
		for (int i=0;i<d2.length;i++)d2[i] = data2[i+offset];
		
		return calculateCorrelation(d1, d2);
		
	}
	
	
	/**
	 * Calculate correlation.
	 * 
	 * @param data1 the first dataset
	 * @param data2 the second dataset
	 * @return the Pearson r-value
	 * @throws SeqMonkException if the two datasets don't have the same number of points in them.
	 */
	public static float calculateCorrelation (float [] data1, float [] data2) {
	
		if (data1.length != data2.length) {
			throw new IllegalArgumentException("Data sets must be the same length when calculating correlation");
		}
		
		float sum12 = 0;
		float sum1 = 0;
		float sum2 = 0;
		float sum1square = 0;
		float sum2square =0;
		
		for (int i=0;i<data1.length;i++) {
			sum12 += data1[i]*data2[i];
			sum1 += data1[i];
			sum2 += data2[i];
			sum1square += data1[i]*data1[i];
			sum2square += data2[i]*data2[i];
		}
		
		float top = sum12 - ((sum1*sum2)/data1.length);
		float bottomRight = sum2square - ((sum2*sum2)/data1.length);
		float bottomLeft = sum1square - ((sum1*sum1)/data1.length);
		float bottom = (float)Math.sqrt(bottomLeft * bottomRight);
		
		
		return top/bottom;
	}
		
}
