/**
 * Copyright 2009- 21 Simon Andrews
 *
 *    This file is part of SeqMonk.
 *
 *    SeqMonk is free software; you can redistribute it and/or modify
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
 *    along with SeqMonk; if not, write to the Free Software
 *    Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */
package uk.ac.babraham.FastQC.Utilities.ImageSaver;

import java.awt.Component;
import java.io.OutputStream;
import java.io.PrintWriter;


/**
 * A utility class which acts as a wrapper for the SVG or PNG generating
 * code which can be used to save (almost) any component which uses the
 * standard Graphics interface to draw itself.
 */
public class SVGImageSaver {

	public static void saveImage (Component c, OutputStream os) {
		PrintWriter pr = new PrintWriter(os);
		SVGGenerator.writeSVG(pr,c);
		pr.flush();
	}
	
	
			
}
