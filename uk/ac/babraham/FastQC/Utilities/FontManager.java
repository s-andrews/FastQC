/**
 * Copyright Copyright 2010-26 Simon Andrews
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
package uk.ac.babraham.FastQC.Utilities;

import java.awt.Font;
import java.awt.FontFormatException;
import java.awt.GraphicsEnvironment;
import java.io.IOException;
import java.io.InputStream;

/**
 * Manages the bundled Liberation Sans font so that FastQC renders
 * consistently without requiring any system-installed fonts.
 */
public class FontManager {

	/** The font family name used in SVG output, with Arial fallback */
	public static final String SVG_FONT_FAMILY = "'Liberation Sans', Arial, Helvetica, sans-serif";

	private static Font regularFont;
	private static Font boldFont;

	private static Font cachedDefaultFont;
	private static Font cachedDefaultBoldFont;

	private static boolean initialized = false;

	/**
	 * Load and register the bundled Liberation Sans font variants.
	 * Safe to call multiple times; only loads once.
	 */
	public static synchronized void initialize() {
		if (initialized) return;

		regularFont = loadFont("/uk/ac/babraham/FastQC/Resources/Fonts/LiberationSans-Regular.ttf");
		boldFont = loadFont("/uk/ac/babraham/FastQC/Resources/Fonts/LiberationSans-Bold.ttf");

		cachedDefaultFont = regularFont.deriveFont(12f);
		cachedDefaultBoldFont = boldFont.deriveFont(12f);

		initialized = true;
	}

	private static Font loadFont(String resourcePath) {
		try (InputStream is = FontManager.class.getResourceAsStream(resourcePath)) {
			if (is == null) {
				System.err.println("Warning: bundled font not found: " + resourcePath);
				return new Font(Font.SANS_SERIF, Font.PLAIN, 12);
			}
			Font font = Font.createFont(Font.TRUETYPE_FONT, is);
			GraphicsEnvironment.getLocalGraphicsEnvironment().registerFont(font);
			return font;
		} catch (FontFormatException | IOException e) {
			System.err.println("Warning: failed to load bundled font: " + resourcePath + " (" + e.getMessage() + ")");
			return new Font(Font.SANS_SERIF, Font.PLAIN, 12);
		}
	}

	/**
	 * Get a font at the given size and style.
	 *
	 * @param style Font.PLAIN or Font.BOLD
	 * @param size  point size
	 * @return the derived font
	 */
	public static Font getFont(int style, float size) {
		initialize();
		Font base = (style == Font.BOLD) ? boldFont : regularFont;
		return base.deriveFont(size);
	}

	/** Default 12pt plain font for graph rendering */
	public static Font defaultFont() {
		initialize();
		return cachedDefaultFont;
	}

	/** Default 12pt bold font for graph rendering */
	public static Font defaultBoldFont() {
		initialize();
		return cachedDefaultBoldFont;
	}
}
