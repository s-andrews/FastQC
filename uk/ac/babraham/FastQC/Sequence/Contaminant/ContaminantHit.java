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
package uk.ac.babraham.FastQC.Sequence.Contaminant;

public class ContaminantHit {

	private Contaminant contaminant;
	private int direction;
	private int length;
	private int percentID;
	
	public static final int FORWARD = 1;
	public static final int REVERSE = 2;
	
	public ContaminantHit (Contaminant contaminant, int direction, int length, int percentID) {
		if (direction == FORWARD || direction == REVERSE) {
			this.direction = direction;
		}
		else {
			throw new IllegalArgumentException("Direction of hit must be FORWARD or REVERSE");
		}
		this.contaminant = contaminant;
		this.length = length;
		this.percentID = percentID;
	}

	
	public Contaminant contaminant () {
		return contaminant;
	}
	
	public int direction () {
		return direction;
	}
	
	public int length () {
		return length;
	}
	
	public int percentID () {
		return percentID;
	}
	
	public String toString () {
		return contaminant.name()+" ("+percentID+"% over "+length+"bp)";
	}
}
