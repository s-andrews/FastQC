/**
 * Copyright Copyright 2014-17 Simon Andrews
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

import java.awt.image.BufferedImage;
import java.io.ByteArrayOutputStream;
import java.io.IOException;
import java.io.OutputStream;

import javax.imageio.ImageIO;

import net.sourceforge.iharder.base64.Base64;

public class ImageToBase64 {

	public static String imageToBase64 (BufferedImage b) {
		ByteArrayOutputStream os = new ByteArrayOutputStream();
		OutputStream b64 = new Base64.OutputStream(os);
		
		try {	
			ImageIO.write(b, "PNG", b64);
		
			return("data:image/png;base64,"+os.toString("UTF-8"));
		}
		catch (IOException e) {
			e.printStackTrace();
			return "Failed";
		}
		
	}

	public static String svgImageToBase64 (String svgdata) {
		

		// We've moved to using the Java.util Base64 encoder which means that
		// SVG output will only work on java v8+ but there was a bug in the 
		// library we were using which caused the last character to get lost
		// some times so this is an easy fix and no one is using java <v8 
		// any more.

		String data = "data:image/svg+xml;base64,"+java.util.Base64.getEncoder().encodeToString(svgdata.getBytes());			

		return(data);
		
	}

	
	
	
}
