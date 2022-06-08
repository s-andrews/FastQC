/**
 * Copyright Copyright 2007-15 Simon Andrews
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

/**
 * A simple class to represent separate values for Red Green and Blue
 * components of a colour.
 */
public class RGB {
	
	public int r;
	public int g;
	public int b;
	
	/**
	 * Instantiates a new RGB colour.
	 * 
	 * @param r RED
	 * @param g GREEN
	 * @param b BLUE
	 */
	public RGB (int r, int g, int b) {
		this.r = r;
		this.g = g;
		this.b = b;
	}
}
