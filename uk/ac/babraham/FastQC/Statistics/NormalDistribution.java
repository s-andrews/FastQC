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

public class NormalDistribution {

	private double mean;
	private double stdev;
	
	public NormalDistribution (double mean, double stdev) {
//		System.out.println("Made distribution with mean "+mean+" and variance "+stdev);
		this.mean = mean;
		this.stdev = stdev;
	}
	
	public double getZScoreForValue (double value) {
		double lhs = 1d/(Math.sqrt(2*Math.PI*stdev*stdev));
		double rhs = Math.pow(Math.E, 0 - (Math.pow(value-mean,2)/(2*stdev*stdev)));
		
		return lhs*rhs;
	}
	
	
//	public static void main (String [] args) {
//		NormalDistribution nd = new NormalDistribution(50, 5);
//		
//		for (int i=0;i<=100;i++) {
//			System.out.println(i+"\t"+nd.getZScoreForValue(i));
//		}
//	}
	
}
