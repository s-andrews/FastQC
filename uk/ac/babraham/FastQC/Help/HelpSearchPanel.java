/**
 * Copyright 2009-15 Simon Andrews
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
package uk.ac.babraham.FastQC.Help;

import java.awt.BorderLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.BorderFactory;
import javax.swing.DefaultListModel;
import javax.swing.JButton;
import javax.swing.JList;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextField;
import javax.swing.ListSelectionModel;
import javax.swing.event.ListSelectionEvent;
import javax.swing.event.ListSelectionListener;

/**
 * The Class HelpSearchPanel.
 */
public class HelpSearchPanel extends JPanel implements ActionListener, ListSelectionListener, Runnable {

	/** The root. */
	private HelpIndexRoot root;
	
	/** The query field. */
	private JTextField queryField;
	
	/** The result list. */
	private JList resultList;
	
	/** The list model. */
	private DefaultListModel listModel;
	
	/** The search button. */
	private JButton searchButton;
	
	/** The dialog. */
	private HelpDialog dialog;
	
	/** The results scroll pane. */
	private JScrollPane resultsScrollPane;
	
	/**
	 * Instantiates a new help search panel.
	 * 
	 * @param root the root
	 * @param dialog the dialog
	 */
	public HelpSearchPanel (HelpIndexRoot root,HelpDialog dialog) {
		this.root = root;
		this.dialog = dialog;
		
		setLayout(new BorderLayout());
		
		JPanel queryPanel = new JPanel();
		queryPanel.setLayout(new BorderLayout());
		queryPanel.setBorder(BorderFactory.createEmptyBorder(2,2,2,2));
		queryField = new JTextField();
		queryField.setActionCommand("search");
		queryField.addActionListener(this);
		queryPanel.add(queryField,BorderLayout.CENTER);
		searchButton = new JButton("Search");
		searchButton.setActionCommand("search");
		searchButton.addActionListener(this);
		queryPanel.add(searchButton,BorderLayout.EAST);
		add(queryPanel,BorderLayout.NORTH);
		
		listModel = new DefaultListModel();
		listModel.addElement("[No search results]");
		resultList = new JList(listModel);
		resultList.addListSelectionListener(this);
		resultList.setSelectionMode(ListSelectionModel.SINGLE_SELECTION);
		resultsScrollPane = new JScrollPane(resultList);
		add(resultsScrollPane,BorderLayout.CENTER);
		
	}

	/* (non-Javadoc)
	 * @see java.awt.event.ActionListener#actionPerformed(java.awt.event.ActionEvent)
	 */
	public void actionPerformed(ActionEvent e) {
		Thread t = new Thread(this);
		t.start();
	}

	/* (non-Javadoc)
	 * @see javax.swing.event.ListSelectionListener#valueChanged(javax.swing.event.ListSelectionEvent)
	 */
	public void valueChanged(ListSelectionEvent lse) {
		Object o = resultList.getSelectedValue();
		if (o != null && o instanceof HelpPage) {
			dialog.DisplayPage((HelpPage)o);
		}
	}

	/* (non-Javadoc)
	 * @see java.lang.Runnable#run()
	 */
	public void run() {
		searchButton.setEnabled(false);
		listModel.removeAllElements();
		if (queryField.getText().trim().length() > 0) {
			HelpPage[] results;
			try {
				results = root.findPagesForTerm(queryField.getText().trim());
			} 
			catch (IOException e) {
				e.printStackTrace();
				searchButton.setEnabled(true);
				return;
			}
			if (results.length > 0) {
				for (int r=0;r<results.length;r++) {
					listModel.addElement(results[r]);
				}
			}
			else {
				listModel.addElement("[No search results]");
			}
		}

		// This stupid rigmarole is because on OSX the updated list
		// just won't show up for some reason.  Removing the list and
		// re-adding it forces it to always show up.
		//
		// It's not even enough to remake the scroll pane.  You have
		// to replace the entire JList.  Aaargh!
		remove(resultsScrollPane);
		revalidate();
		resultList = new JList(listModel);
		resultList.addListSelectionListener(this);
		resultsScrollPane = new JScrollPane(resultList);
		add(resultsScrollPane,BorderLayout.CENTER);
		revalidate();
		repaint();
		
		searchButton.setEnabled(true);
	}
}
