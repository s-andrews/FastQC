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

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;

import javax.swing.JLabel;
import javax.swing.JPanel;

@SuppressWarnings("serial")
public class WelcomePanel extends JPanel {
	private static final long serialVersionUID = 1L;

	public WelcomePanel () {
		setLayout(new GridBagLayout());
		GridBagConstraints gbc = new GridBagConstraints();
		
		gbc.gridx=1;
		gbc.gridy=1;
		gbc.weightx=0.5;
		gbc.weighty=0.99;
		
		add(new JPanel(),gbc);
		gbc.gridy++;
		gbc.weighty=0.01;
		
		gbc.insets = new Insets(10, 10, 10, 10);
		gbc.fill = GridBagConstraints.NONE;
		
		add(new FastQCTitlePanel(),gbc);
		
		gbc.gridy++;
		gbc.weighty=0.5;
		
		add(new JLabel("Use File > Open to select the sequence file you want to check"),gbc);
		
		gbc.gridy++;
		gbc.weighty=0.99;
		add(new JPanel(),gbc);
		
	}
}
