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

import javax.swing.*;

import uk.ac.babraham.FastQC.FastQCApplication;

import java.awt.*;
import java.awt.event.*;

/**
 * Shows the generic about dialog giving details of the current version
 * and copyright assignments.  This is just a thin shell around the 
 * SeqMonkTitlePanel which actually holds the relevant information and
 * which is also used on the welcome screen.
 */
public class AboutDialog extends JDialog {
	private static final long serialVersionUID = 1L;

    /**
     * Instantiates a new about dialog.
     * 
     * @param a The SeqMonk application.
     */
    public AboutDialog(FastQCApplication a) {
    	super(a);
        setTitle("About FastQC...");  
        Container cont = getContentPane();
        cont.setLayout(new BorderLayout());
        
        add(new FastQCTitlePanel(),BorderLayout.CENTER);

        JPanel buttonPanel = new JPanel();
        
        JButton closeButton = new JButton("Close");
        getRootPane().setDefaultButton(closeButton);
        closeButton.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent arg0) {
                setVisible(false);
                dispose();
            }
        });
        buttonPanel.add(closeButton);
        
        cont.add(buttonPanel,BorderLayout.SOUTH);
        
        setSize(650,230);
        setLocationRelativeTo(a);
        setResizable(false);
        setVisible(true);
    }
    
}
