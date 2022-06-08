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
package uk.ac.babraham.FastQC.Dialogs;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Font;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.awt.RenderingHints;

import javax.swing.BorderFactory;
import javax.swing.ImageIcon;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JTextField;

import uk.ac.babraham.FastQC.FastQCApplication;

/**
 * The Class SeqMonkTitlePanel.
 */

public class FastQCTitlePanel extends JPanel {
	private static final long serialVersionUID = 1L;

	/**
	 * Provides a small panel which gives details of the FastQC version
	 * and copyright.  Used in both the welcome panel and the about dialog.
	 */
	public FastQCTitlePanel () {
		setLayout(new BorderLayout(5,1));

		ImageIcon logo = new ImageIcon(ClassLoader.getSystemResource("uk/ac/babraham/FastQC/Resources/fastqc_icon_100.png"));
		JPanel logoPanel = new JPanel();
		logoPanel.add(new JLabel("",logo,JLabel.CENTER));
		logoPanel.setBorder(BorderFactory.createEmptyBorder(3, 3, 3, 3));
		add(logoPanel,BorderLayout.WEST);
		JPanel c = new JPanel();
		c.setLayout(new GridBagLayout());

		GridBagConstraints constraints = new GridBagConstraints();
		constraints.gridx=1;
		constraints.gridy=1;
		constraints.weightx = 1;
		constraints.weighty=1;
		constraints.insets = new Insets(3, 3, 0, 0);
		constraints.fill = GridBagConstraints.NONE;

		JLabel program = new SmoothJLabel("FastQC High Throughput Sequence QC Report",JLabel.CENTER);
		program.setFont(new Font("Dialog",Font.BOLD,18));
		program.setForeground(new Color(200,0,0));
		c.add(program,constraints);

		constraints.gridy++;
		JLabel version = new SmoothJLabel("Version: "+FastQCApplication.VERSION, JLabel.CENTER);
		version.setFont(new Font("Dialog",Font.BOLD,15));
		version.setForeground(new Color(0,0,200));
		c.add(version,constraints);

		constraints.gridy++;
		// Use a text field so they can copy this
		JTextField website = new JTextField(" www.bioinformatics.babraham.ac.uk/projects/ ");
		website.setFont(new Font("Dialog",Font.PLAIN,14));
		website.setEditable(false);
		website.setBorder(null);
		website.setOpaque(false);
		website.setHorizontalAlignment(JTextField.CENTER);
		c.add(website,constraints);
		constraints.gridy++;

		JLabel copyright = new JLabel("\u00a9 Simon Andrews, Pierre Lindenbaum, Brian Howard, Phil Ewels 2011-22,", JLabel.CENTER);
		copyright.setFont(new Font("Dialog",Font.PLAIN,14));
		c.add(copyright,constraints);
		constraints.gridy++;
		
		JLabel copyright2 = new JLabel("HTSJDK BAM/SAM reader v2.24.1 from Samtools 2022", JLabel.CENTER);
		copyright2.setFont(new Font("Dialog",Font.PLAIN,10));
		c.add(copyright2,constraints);
		constraints.gridy++;
		
		JLabel copyright3 = new JLabel("BZip decompression \u00a9Matthew J. Francis, 2011", JLabel.CENTER);
		copyright3.setFont(new Font("Dialog",Font.PLAIN,10));
		c.add(copyright3,constraints);
		constraints.gridy++;
		
		JLabel copyright4 = new JLabel("Base64 encoding \u00a9Robert Harder, 2012", JLabel.CENTER);
		copyright4.setFont(new Font("Dialog",Font.PLAIN,10));
		c.add(copyright4,constraints);
		constraints.gridy++;
		
		JLabel copyright5 = new JLabel("Java HDF5 reader \u00a9ETH, CISD and SIS, 2007-14", JLabel.CENTER);
		copyright5.setFont(new Font("Dialog",Font.PLAIN,10));
		c.add(copyright5,constraints);

		add(c,BorderLayout.CENTER);
	}
	
	/**
	 * A JLabel with anti-aliasing enabled.  Takes the same constructor
	 * arguments as JLabel
	 */
	private class SmoothJLabel extends JLabel {
		
		/**
		 * Creates a new label
		 * 
		 * @param text The text
		 * @param position The JLabel constant position for alignment
		 */
		public SmoothJLabel (String text, int position) {
			super(text,position);
		}
		
		/* (non-Javadoc)
		 * @see javax.swing.JComponent#paintComponent(java.awt.Graphics)
		 */
		public void paintComponent (Graphics g) {
			if (g instanceof Graphics2D) {
				((Graphics2D)g).setRenderingHint(RenderingHints.KEY_ANTIALIASING, RenderingHints.VALUE_ANTIALIAS_ON);
			}
			super.paintComponent(g);
		}

	}
	
}
