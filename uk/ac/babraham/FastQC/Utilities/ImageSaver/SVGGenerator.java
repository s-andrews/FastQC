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

import java.awt.Color;
import java.awt.Component;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Image;
import java.awt.Rectangle;
import java.awt.Shape;
import java.awt.image.BufferedImage;
import java.awt.image.ImageObserver;
import java.text.AttributedCharacterIterator;

import javax.swing.RepaintManager;

/**
 * SVGGenerator is an implementation of the Graphics interface which can be
 * passed to any paint method and will generate an SVG version of the display.
 * 
 * This is not a complete implementation of the graphics interface - just the
 * primitive drawing parts.  It doesn't handle drawing bitmaps and it doesn't
 * implement any of the Graphics2D extensions so you should test it against
 * any class you want to convert into SVG.
 * 
 * Generally this class shouldn't be called directly.  You should use the
 * ImageSaver class which allows the selection of the file to save and the
 * type of image to generate.  This class will be called by ImageSaver.
 */
public class SVGGenerator {

	/**
	 * This class is used to generate an SVG version of displays
	 * created within the program.  Although it extends Graphics
	 * so we can pass it to a paint routine it's not trying to
	 * reimplement the whole of the Graphics interface, just the
	 * bits we actually use.	
	 */
		
	private StringBuffer sb;
	
	private Graphics g;
	private int cWidth = 0;
	private int cHeight = 0;
	
	/*
	 * These are the methods we know we actually use so we need
	 * to implement them in our generator
	 */
	
	/**
	 * This is the method used to create an SVG representation
	 * of a component.  The returned string contains valid SVG
	 * 
	 * @param c The component to convert
	 * @return An SVG representation of the component
	 */
	public static String writeSVG (Component c) {
		
		StringBuffer sb = new StringBuffer();
		
		/*
		 * Before using our Graphics class we need to disable double
		 * buffering on the component.  If we don't do this then we
		 * just get an image from the offscreen buffer to draw into
		 * our Graphics object - we never see the individual method
		 * calls to the Graphics interface.
		 * 
		 * This only affects Windows/Linux where java does the buffering
		 * internally.  OSX does it via the window manager so this works
		 * without changing anything.
		 */
				
		boolean doubleBuffered = RepaintManager.currentManager(c).isDoubleBufferingEnabled();
		if (doubleBuffered) {
			RepaintManager.currentManager(c).setDoubleBufferingEnabled(false);
		}
		
		new SVGGenerator(sb,c);
		
		if (doubleBuffered) {
			RepaintManager.currentManager(c).setDoubleBufferingEnabled(true);
		}
		
		return(sb.toString());

	}
	
	/**
	 * Instantiates a new generator.  Not used externally - all external calls
	 * to this class should go via the static convert to SVG method.
	 * 
	 * @param c The component to convert
	 */
	private SVGGenerator (StringBuffer sb, Component c) {
			
		this.sb = sb;
		
		cWidth = c.getWidth();
		cHeight = c.getHeight();
		
		// Some methods would be a pain to reimplement ourselves and
		// aren't used directly by us (eg providing font metrics).  We
		// therefore create a graphics object from a BufferedImage and
		// use that to do any difficult stuff we can't be bothered with.
		
		BufferedImage b = new BufferedImage(c.getWidth(),c.getHeight(),BufferedImage.TYPE_INT_RGB);
		g = b.getGraphics();			

		
		sb.append("<?xml version=\"1.0\" standalone=\"no\"?>\n");
		sb.append("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" \"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">\n");
		sb.append("<svg width=\"");
		sb.append(c.getWidth());
		sb.append("\" height=\"");
		sb.append(c.getHeight());
		sb.append("\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">\n");

		c.paint(new SVGGraphics());

		sb.append("</svg>\n");
		
	}
	
	
	
		
	/**
	 * The Actual SVG implementation of Graphics
	 */
	private class SVGGraphics extends Graphics {
		
		/**
		 * SVGGraphics is an implementation of the Graphics
		 * abstract class.  When used it is required to copy
		 * itself and inherit translation and clipping properties
		 * from its parent.  We don't worry about clipping since
		 * we're only going to get drawn once.  We also store font
		 * and color although I'm not sure if that's necessary.
		 */
		
		/** The current offset on the x-axis **/
		private int translateX = 0;		

		/** The current offset on the y-axis **/
		private int translateY = 0;
		
		private Color color = Color.BLACK;
		
		/** The font. */
		private Font font = new Font("Default",Font.PLAIN,12);
		

		/**
		 * Instantiates a new SVG graphics object.
		 */
		public SVGGraphics () {
		}
		
		/**
		 * SVG Graphics has to recursively call itself.  This shouldn't be
		 * called externally.
		 * 
		 * @param parent The calling SVGGraphics instance
		 */
		public SVGGraphics(SVGGraphics parent) {
			this.translateX = parent.translateX;
			this.translateY = parent.translateY;
			this.color = parent.color;
			this.font = parent.font;
		}
	
		/**
		 * Corrects an x value with the current translation value
		 * 
		 * @param x The raw x value
		 * @return The corrected x value
		 */
		private int correctX (int x){
			int newX = x + translateX;
			if (newX<0) newX=0;
			if (newX>cWidth) newX=cWidth;
			
			return newX;
		}
	
		/**
		 * Corrects a y value with the current translation value
		 * 
		 * @param y The raw y value
		 * @return The corrected y value
		 */
		private int correctY (int y) {
			int newY = y + translateY;
			if (newY<0) newY=0;
			if (newY>cHeight) newY=cHeight;

			return newY;
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawLine(int, int, int, int)
		 */
		@Override
		public void drawLine(int x1, int y1, int x2, int y2) {
			x1 = correctX(x1);
			y1 = correctY(y1);
			x2 = correctX(x2);
			y2 = correctY(y2);
			
			sb.append("<line x1=\"");
			sb.append(x1);
			sb.append("\" y1=\"");
			sb.append(y1);
			sb.append("\" x2=\"");
			sb.append(x2);
			sb.append("\" y2=\"");
			sb.append(y2);
			sb.append("\" stroke=\"rgb(");
			appendColor();
			sb.append(")\" stroke-width=\"1\"/>\n");
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawOval(int, int, int, int)
		 */
		@Override
		public void drawOval(int x, int y, int width, int height) {
			
			x = correctX(x);
			y = correctY(y);
			
			if (x+width > cWidth) {
				width = cWidth-x;
			}
	
			if (y+height > cHeight) {
				height = cHeight-y;
			}
			
			// SVG does ovals from the middle of the circle rather
			// than the edge so we need to do some conversion.
			
			x+= width/2;
			y+=height/2;
			width/=2;
			height/=2;
			
			sb.append("<ellipse cx=\"");
			sb.append(x);
			sb.append("\" cy=\"");
			sb.append(y);
			sb.append("\" rx=\"");
			sb.append(width);
			sb.append("\" ry=\"");
			sb.append(height);
			sb.append("\" style=\"fill:none;stroke:rgb(");
			appendColor();
			sb.append(");stroke-width:1\"/>\n");
		}
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#fillOval(int, int, int, int)
		 */
		@Override
		public void fillOval(int x, int y, int width, int height) {
			
			x = correctX(x);
			y = correctY(y);
			
			if (x+width > cWidth) {
				width = cWidth-x;
			}
	
			if (y+height > cHeight) {
				height = cHeight-y;
			}
			
			// SVG does ovals from the middle of the circle rather
			// than the edge so we need to do some conversion.
			
			x+= width/2;
			y+=height/2;
			width/=2;
			height/=2;
			
			sb.append("<ellipse cx=\"");
			sb.append(x);
			sb.append("\" cy=\"");
			sb.append(y);
			sb.append("\" rx=\"");
			sb.append(width);
			sb.append("\" ry=\"");
			sb.append(height);
			sb.append("\" style=\"fill:rgb(");
			appendColor();
			sb.append(");stroke:none\"/>\n");
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawString(java.lang.String, int, int)
		 */
		@Override
		public void drawString(String string, int x, int y) {
			
			x = correctX(x);
			y = correctY(y);
			
			if (string.indexOf("&")>=0) {
				string = "<![CDATA[" + string.replaceAll("]]>", "]]>]]><![CDATA[")+ "]]>";
			}
			else {
				string = string.replaceAll(">", "&gt;");
				string = string.replaceAll("<", "&lt;");
			}
			
			// This doesn't yet take the font family into account
			// as there isn't always an easy way to map between
			// java names and SVG names.  Will look at this again
			// later.
			
			sb.append("<text x=\"");
			sb.append(x);
			sb.append("\" y=\"");
			sb.append(y);
			sb.append("\" fill=\"rgb(");
			appendColor();
			sb.append(")\" font-family=\"Arial\" font-size=\"");
			sb.append(font.getSize());
			sb.append("\">");
			sb.append(string);
			sb.append("</text>\n");
		}
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#setColor(java.awt.Color)
		 */
		public void setColor(Color color) {
			this.color = color;
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#setFont(java.awt.Font)
		 */
		public void setFont(Font font) {
			this.font = font;
		}
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawRect(int, int, int, int)
		 */
		public void drawRect(int x, int y, int width, int height) {
			drawRoundRect(x, y, width, height, 0, 0);
		}
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawRoundRect(int, int, int, int, int, int)
		 */
		public void drawRoundRect(int x, int y, int width, int height, int arcWidth,int arcHeight) {
	
			x = correctX(x);
			y = correctY(y);
			
			if (x+width > cWidth) {
				width = cWidth-x;
			}
	
			if (y+height > cHeight) {
				height = cHeight-y;
			}
			
			sb.append("<rect width=\"");
			sb.append(width);
			sb.append("\" height=\"");
			sb.append(height);
			sb.append("\"");
			sb.append(" x=\"");
			sb.append(x);
			sb.append("\" y=\"");
			sb.append(y);
			sb.append("\" rx=\"");
			sb.append(arcWidth);
			sb.append("\" ry=\"");
			sb.append(arcHeight);
			sb.append("\" style=\"fill:none");
			sb.append(";stroke-width:1;stroke:rgb(");
					
			appendColor();
					
			sb.append(")\"/>\n");
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#fillRect(int, int, int, int)
		 */
		public void fillRect(int x, int y, int width, int height) {
			fillRoundRect(x, y, width, height, 0, 0);	
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#fillRoundRect(int, int, int, int, int, int)
		 */
		public void fillRoundRect(int x, int y, int width, int height, int arcWidth,int arcHeight) {
	
			x = correctX(x);
			y = correctY(y);
						
			if (x+width > cWidth) {
				width = cWidth-x;
			}
	
			if (y+height > cHeight) {
				height = cHeight-y;
			}
			
			sb.append("<rect width=\"");
			sb.append(width);
			sb.append("\" height=\"");
			sb.append(height);
			sb.append("\"");
			sb.append(" x=\"");
			sb.append(x);
			sb.append("\" y=\"");
			sb.append(y);
			if (arcWidth > 0 || arcHeight > 0) {
				sb.append("\" rx=\"");
				sb.append(arcWidth);
				sb.append("\" ry=\"");
				sb.append(arcHeight);
			}
			sb.append("\" style=\"fill:rgb(");
			
			appendColor();
	
			sb.append(");stroke:none\"/>\n");
		}
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#toString()
		 */
		public String toString () {
			return "SVGGraphics";
		}
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#translate(int, int)
		 */
		public void translate(int translateX, int translateY) {
			/**
			 * It's important to note that new translations are
			 * layered on top of previous translations and do
			 * not replace them.  We got this wrong previously
			 * and it made a right mess!
			 */
			
			this.translateX += translateX;
			this.translateY += translateY;
		}
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#create()
		 */
		public Graphics create() {
			return new SVGGraphics(this);
		}
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#getFont()
		 */
		public Font getFont() {
			return font;
		}
		
		/**
		 * Append color.
		 */
		private void appendColor () {
			sb.append(color.getRed());
			sb.append(",");
			sb.append(color.getGreen());
			sb.append(",");
			sb.append(color.getBlue());
		}
	
		/*
		 * The following methods are required by the graphics interface
		 * but we don't actually use them so we're just leaving them as
		 * stubs.
		 */
		
		/* (non-Javadoc)
		 * @see java.awt.Graphics#clearRect(int, int, int, int)
		 */
		public void clearRect(int arg0, int arg1, int arg2, int arg3) {}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#clipRect(int, int, int, int)
		 */
		public void clipRect(int x, int y, int width, int height) {}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#copyArea(int, int, int, int, int, int)
		 */
		public void copyArea(int arg0, int arg1, int arg2, int arg3, int arg4,int arg5) {}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#dispose()
		 */
		public void dispose() {}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawArc(int, int, int, int, int, int)
		 */
		public void drawArc(int arg0, int arg1, int arg2, int arg3, int arg4,int arg5) {}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawImage(java.awt.Image, int, int, java.awt.image.ImageObserver)
		 */
		public boolean drawImage(Image img, int x, int y, ImageObserver obs) {
			System.out.println("Height="+img.getHeight(obs)+" Width="+img.getWidth(obs));
			return g.drawImage(img, x, y, obs);
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawImage(java.awt.Image, int, int, java.awt.Color, java.awt.image.ImageObserver)
		 */
		public boolean drawImage(Image arg0, int arg1, int arg2, Color arg3,ImageObserver arg4) {
			System.out.println("Draw image 2");
			return true;
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawImage(java.awt.Image, int, int, int, int, java.awt.image.ImageObserver)
		 */
		public boolean drawImage(Image arg0, int arg1, int arg2, int arg3,int arg4, ImageObserver arg5) {
			System.out.println("Draw image 3");
			return true;
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawImage(java.awt.Image, int, int, int, int, java.awt.Color, java.awt.image.ImageObserver)
		 */
		public boolean drawImage(Image arg0, int arg1, int arg2, int arg3,int arg4, Color arg5, ImageObserver arg6) {
			System.out.println("Draw image 4");
			return true;
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawImage(java.awt.Image, int, int, int, int, int, int, int, int, java.awt.image.ImageObserver)
		 */
		public boolean drawImage(Image arg0, int arg1, int arg2, int arg3,int arg4, int arg5, int arg6, int arg7, int arg8, ImageObserver arg9) {
			System.out.println("Draw image 5");
			return true;
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawImage(java.awt.Image, int, int, int, int, int, int, int, int, java.awt.Color, java.awt.image.ImageObserver)
		 */
		public boolean drawImage(Image arg0, int arg1, int arg2, int arg3,int arg4, int arg5, int arg6, int arg7, int arg8, Color arg9,ImageObserver arg10) {
			System.out.println("Draw image 6");
			return true;
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawPolygon(int[], int[], int)
		 */
		public void drawPolygon(int[] xPoints, int[] yPoints, int pointsToUse) {
			
			int [] correctedXPoints = new int[pointsToUse];
			int [] correctedYPoints = new int[pointsToUse];
			
			for (int i=0;i<pointsToUse;i++) {
				correctedXPoints[i] = correctX(xPoints[i]);
				correctedYPoints[i] = correctY(yPoints[i]);
			}
						
			sb.append("<polygon points=\"");

			for (int i=0;i<pointsToUse;i++) {
				if (i != 0) {
					sb.append(",");
				}
				sb.append(correctedXPoints[i]);
				sb.append(",");
				sb.append(correctedYPoints[i]);
			}
			
			sb.append("\" style=\"stroke-width:1;stroke:rgb(");
			
			appendColor();
	
			sb.append(");fill:none\"/>\n");
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawPolyline(int[], int[], int)
		 */
		public void drawPolyline(int[] arg0, int[] arg1, int arg2) {
			System.out.println("Draw polygon 2");
		}
	
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#drawString(java.text.AttributedCharacterIterator, int, int)
		 */
		public void drawString(AttributedCharacterIterator arg0, int arg1, int arg2) {
			System.out.println("Draw complex string");
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#fillArc(int, int, int, int, int, int)
		 */
		public void fillArc(int arg0, int arg1, int arg2, int arg3, int arg4,int arg5) {
			System.out.println("Fill arc");
		}
	
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#fillPolygon(int[], int[], int)
		 */
		public void fillPolygon(int[] xPoints, int[] yPoints, int pointsToUse) {
			
			int [] correctedXPoints = new int[pointsToUse];
			int [] correctedYPoints = new int[pointsToUse];
			
			for (int i=0;i<pointsToUse;i++) {
				correctedXPoints[i] = correctX(xPoints[i]);
				correctedYPoints[i] = correctY(yPoints[i]);
			}
						
			sb.append("<polygon points=\"");

			for (int i=0;i<pointsToUse;i++) {
				if (i != 0) {
					sb.append(",");
				}
				sb.append(correctedXPoints[i]);
				sb.append(",");
				sb.append(correctedYPoints[i]);
			}
			
			sb.append("\" style=\"fill:rgb(");
			
			appendColor();
	
			sb.append(");stroke:none\"/>\n");

		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#getClip()
		 */
		public Shape getClip() {
			return g.getClip();
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#getClipBounds()
		 */
		public Rectangle getClipBounds() {
			return g.getClipBounds();
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#getColor()
		 */
		public Color getColor() {
			return color;
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#getFontMetrics(java.awt.Font)
		 */
		public FontMetrics getFontMetrics(Font f) {
			return g.getFontMetrics(f);
		}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#setClip(java.awt.Shape)
		 */
		public void setClip(Shape s) {}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#setClip(int, int, int, int)
		 */
		public void setClip(int x, int y, int width, int height) {}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#setPaintMode()
		 */
		public void setPaintMode() {}
	
		/* (non-Javadoc)
		 * @see java.awt.Graphics#setXORMode(java.awt.Color)
		 */
		public void setXORMode(Color arg0) {}
	}
	
}
