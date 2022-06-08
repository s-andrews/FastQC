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
package uk.ac.babraham.FastQC;

import java.awt.BorderLayout;
import java.io.File;
import java.io.IOException;
import java.util.Vector;

import javax.swing.ImageIcon;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.UIManager;
import javax.swing.filechooser.FileFilter;

import uk.ac.babraham.FastQC.Analysis.AnalysisRunner;
import uk.ac.babraham.FastQC.Analysis.OfflineRunner;
import uk.ac.babraham.FastQC.Dialogs.WelcomePanel;
import uk.ac.babraham.FastQC.FileFilters.BAMFileFilter;
import uk.ac.babraham.FastQC.FileFilters.CasavaFastQFileFilter;
import uk.ac.babraham.FastQC.FileFilters.FastQFileFilter;
import uk.ac.babraham.FastQC.FileFilters.MappedBAMFileFilter;
import uk.ac.babraham.FastQC.FileFilters.SequenceFileFilter;
import uk.ac.babraham.FastQC.Modules.ModuleFactory;
import uk.ac.babraham.FastQC.Modules.QCModule;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Results.ResultsPanel;
import uk.ac.babraham.FastQC.Sequence.SequenceFactory;
import uk.ac.babraham.FastQC.Sequence.SequenceFile;
import uk.ac.babraham.FastQC.Sequence.SequenceFormatException;
import uk.ac.babraham.FastQC.Utilities.CasavaBasename;
import uk.ac.babraham.FastQC.Utilities.NanoporeBasename;

public class FastQCApplication extends JFrame {	
	
	public static final String VERSION = "0.11.10.devel";
	
	private JTabbedPane fileTabs;
	private WelcomePanel welcomePanel;
	private File lastUsedDir = null;
	
	public FastQCApplication () {
			setTitle("FastQC");
			setIconImage(new ImageIcon(ClassLoader.getSystemResource("uk/ac/babraham/FastQC/Resources/fastqc_icon.png")).getImage());
			setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
	//		setSize(1280, 720);
			setSize(800,600);
			setLocationRelativeTo(null);
			
			welcomePanel = new WelcomePanel();
			
			fileTabs = new JTabbedPane(JTabbedPane.TOP);
			setContentPane(welcomePanel);
			
			setJMenuBar(new FastQCMenuBar(this));
			
		}

	public void close () {
		if (fileTabs.getSelectedIndex() >=0) {
			fileTabs.remove(fileTabs.getSelectedIndex());
		}
		if (fileTabs.getTabCount() == 0) {
			setContentPane(welcomePanel);
			validate();
			repaint();
		}
	}
	
	public void closeAll () {
		fileTabs.removeAll();
		setContentPane(welcomePanel);
		validate();
		repaint();
	}
	
	public void openFile () {
		JFileChooser chooser;
		
		if (lastUsedDir == null) {
			chooser = new JFileChooser();
		}
		else {
			chooser = new JFileChooser(lastUsedDir);
		}
		chooser.setMultiSelectionEnabled(true);
		SequenceFileFilter sff = new SequenceFileFilter();
		chooser.addChoosableFileFilter(sff);
		chooser.addChoosableFileFilter(new FastQFileFilter());
		chooser.addChoosableFileFilter(new CasavaFastQFileFilter());
		chooser.addChoosableFileFilter(new BAMFileFilter());
		chooser.addChoosableFileFilter(new MappedBAMFileFilter());
		chooser.setFileFilter(sff);
		int result = chooser.showOpenDialog(this);
		if (result == JFileChooser.CANCEL_OPTION) return;
	
		// See if they forced a file format
		FileFilter chosenFilter = chooser.getFileFilter();
		if (chosenFilter instanceof FastQFileFilter) {
			FastQCConfig.getInstance().setSequenceFormat("fastq");
		}
		if (chosenFilter instanceof CasavaFastQFileFilter) {
			FastQCConfig.getInstance().setSequenceFormat("fastq");
			FastQCConfig.getInstance().setCasavaMode(true);

		}
		else if (chosenFilter instanceof BAMFileFilter) {
			FastQCConfig.getInstance().setSequenceFormat("bam");
		}
		else if (chosenFilter instanceof MappedBAMFileFilter) {
			FastQCConfig.getInstance().setSequenceFormat("bam_mapped");
			System.setProperty("fastqc.sequence_format", "bam_mapped");
		}
		
		// If we're still showing the welcome panel switch this out for
		// the file tabs panel
		if (fileTabs.getTabCount() == 0) {
			setContentPane(fileTabs);
			validate();
			repaint();
		}
		
		File [] files = chooser.getSelectedFiles();
		
		if (FastQCConfig.getInstance().nano) {
			// Some of the files we've been passed might be directories and we would need to
			// list the fast5 files within those directories.
			
			Vector<File> keptFiles = new Vector<File>();
			
			for (int f=0;f<files.length;f++) {
				if (files[f].isDirectory()) {
					File [] fast5files = files[f].listFiles();
					for (int i=0;i<fast5files.length;i++) {
						if (fast5files[i].getName().endsWith(".fast5")) {
							keptFiles.add(fast5files[i]);
						}
					}
				}
				else {
					keptFiles.add(files[f]);
				}
			}
			
			files = keptFiles.toArray(new File[0]);
			
		}
		
		
		File [][] fileGroups;
		
		// See if we need to group together files from a casava group
		if (FastQCConfig.getInstance().casava) {
			fileGroups = CasavaBasename.getCasavaGroups(files);
		}
		else if (FastQCConfig.getInstance().nano) {
			fileGroups = NanoporeBasename.getNanoporeGroups(files);
		}
		else {
			fileGroups = new File [files.length][1];
			for (int f=0;f<files.length;f++) {
				fileGroups[f][0] = files[f];
			}
		}

	
		for (int i=0;i<fileGroups.length;i++) {
			File [] filesToProcess = fileGroups[i];
			lastUsedDir = filesToProcess[0].getParentFile();
			SequenceFile sequenceFile;
			
			
			try {
				sequenceFile = SequenceFactory.getSequenceFile(filesToProcess);
			}
			catch (SequenceFormatException e) {
				JPanel errorPanel = new JPanel();
				errorPanel.setLayout(new BorderLayout());
				errorPanel.add(new JLabel("File format error: "+e.getLocalizedMessage(), JLabel.CENTER),BorderLayout.CENTER);
				fileTabs.addTab(filesToProcess[0].getName(), errorPanel);
				e.printStackTrace();
				continue;
			}
			catch (IOException e) {
				System.err.println("File broken");
				e.printStackTrace();
				JOptionPane.showMessageDialog(this, "Couldn't read file:"+e.getLocalizedMessage(), "Error reading file", JOptionPane.ERROR_MESSAGE);
				continue;
			}
					
			AnalysisRunner runner = new AnalysisRunner(sequenceFile);
			ResultsPanel rp = new ResultsPanel(sequenceFile);
			runner.addAnalysisListener(rp);
			fileTabs.addTab(sequenceFile.name(), rp);
			

			QCModule [] module_list = ModuleFactory.getStandardModuleList();
	
			runner.startAnalysis(module_list);
		}
	}

	public void saveReport () {
		JFileChooser chooser;
		
		if (lastUsedDir == null) {
			chooser = new JFileChooser();
		}
		else {
			chooser = new JFileChooser(lastUsedDir);
		}
		
		if (fileTabs.getSelectedComponent() == null) {
			JOptionPane.showMessageDialog(this, "No FastQ files are open yet", "Can't save report", JOptionPane.ERROR_MESSAGE);
			return;
		}
		chooser.setSelectedFile(new File(((ResultsPanel)fileTabs.getSelectedComponent()).sequenceFile().getFile().getName().replaceAll("stdin:","").replaceAll(".gz$","").replaceAll(".bz2$","").replaceAll(".txt$","").replaceAll(".fastq$", "").replaceAll(".fq$", "").replaceAll(".sam$", "").replaceAll(".bam$", "")+"_fastqc.html"));
		chooser.setMultiSelectionEnabled(false);
		chooser.setFileFilter(new FileFilter() {
		
			public String getDescription() {
				return "HTML files";
			}
		
			public boolean accept(File f) {
				if (f.isDirectory() || f.getName().toLowerCase().endsWith(".html")) {
					return true;
				}
				else {
					return false;
				}
			}
		
		});
	
		File reportFile;
		while (true) {
			int result = chooser.showSaveDialog(this);
			if (result == JFileChooser.CANCEL_OPTION) return;
			
			reportFile = chooser.getSelectedFile();
			if (! reportFile.getName().toLowerCase().endsWith(".html")) {
				reportFile = new File(reportFile.getAbsoluteFile()+".html");
			}
			
			// Check if we're overwriting something
			if (reportFile.exists()) {
				int reply = JOptionPane.showConfirmDialog(this, reportFile.getName()+" already exists.  Overwrite?", "Overwrite existing file?", JOptionPane.YES_NO_OPTION);
				if (reply == JOptionPane.NO_OPTION) {
					continue;
				}
				else {
					break;
				}
			}
			else {
				break;
			}
		}
		
		
		ResultsPanel selectedPanel = (ResultsPanel)fileTabs.getSelectedComponent();
		
		try {
			new HTMLReportArchive(selectedPanel.sequenceFile(), selectedPanel.modules(), reportFile);
		} 
		catch (Exception e) {
			JOptionPane.showMessageDialog(this, "Failed to create archive: "+e.getMessage(), "Error", JOptionPane.ERROR_MESSAGE);
			e.printStackTrace();
		}
	}

	public static void main(String[] args) {
		
		// See if we just have to print out the version
		if (System.getProperty("fastqc.show_version") != null && System.getProperty("fastqc.show_version").equals("true")) {
			System.out.println("FastQC v"+VERSION);
			System.exit(0);
		}
		
		if (args.length > 0) {
			// Set headless to true so we don't get problems
			// with people working without an X display.
			System.setProperty("java.awt.headless", "true");
			
			// We used to default to unzipping the zip file in 
			// non-interactive runs.  As we now save an HTML
			// report at the top level we no longer do this
			// so unzip is false unless explicitly set to be true.
						
			if (FastQCConfig.getInstance().do_unzip == null) {
				FastQCConfig.getInstance().do_unzip = false;
			}
			
			new OfflineRunner(args);
			System.exit(0);
		}
		
		else {
			// Recent java themes for linux are just horribly broken with missing
			// bits of UI.  We're therefore not going to set a native look if
			// we're on linux.  See seqmonk bug #95 for details.
			
			try {
				if (! System.getProperty("os.name").toLowerCase().contains("linux")) {
					UIManager.setLookAndFeel(UIManager.getSystemLookAndFeelClassName());
				}
			} catch (Exception e) {}
			
	
			// The interactive default is to not uncompress the
			// reports after they have been generated
			if (FastQCConfig.getInstance().do_unzip == null) {
				FastQCConfig.getInstance().do_unzip = false;
			}
	
			FastQCApplication app = new FastQCApplication();
	
			app.setVisible(true);
		}
	}	

}
