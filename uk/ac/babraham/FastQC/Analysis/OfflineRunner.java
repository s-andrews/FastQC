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
package uk.ac.babraham.FastQC.Analysis;

import java.io.File;
import java.io.IOException;
import java.util.Vector;
import java.util.concurrent.atomic.AtomicInteger;

import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Modules.ModuleFactory;
import uk.ac.babraham.FastQC.Modules.QCModule;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Sequence.SequenceFactory;
import uk.ac.babraham.FastQC.Sequence.SequenceFile;
import uk.ac.babraham.FastQC.Utilities.CasavaBasename;
import uk.ac.babraham.FastQC.Utilities.NanoporeBasename;

public class OfflineRunner implements AnalysisListener {
	
	private AtomicInteger filesRemaining;
	private boolean showUpdates = true;
	
	// We'll set a flag which will tell us if any of the sequences failed
	// so we can exit with an error state so the calling process can tell
	// that something went wrong.
	
	private boolean somethingFailed = false;
	
	public OfflineRunner (String [] filenames) {	
		
		// See if we need to show updates
		showUpdates = !FastQCConfig.getInstance().quiet;
		
		Vector<File> files = new Vector<File>();
		
		// We make a special case if they supply a single filename
		// which is stdin.  In this case we'll take data piped to us
		// rather than trying to read the actual file.  We'll also
		// skip the existence check.
				
		if (filenames.length == 1 && filenames[0].startsWith("stdin")) {
			files.add(new File(filenames[0]));
		}
		else {
			for (int f=0;f<filenames.length;f++) {
				File file = new File(filenames[f]);

				if (!file.exists() || ! file.canRead()) {
					System.err.println("Skipping '"+filenames[f]+"' which didn't exist, or couldn't be read");
					continue;
				}

				if (FastQCConfig.getInstance().nano && file.isDirectory()) {
					File [] fast5files = file.listFiles();
					for (int i=0;i<fast5files.length;i++) {
						if (fast5files[i].getName().endsWith(".fast5")) {
							files.add(fast5files[i]);
						}
					}
					
					// In newer nanopore software instances the fast5 files are
					// put into subdirectories of the main one specified so we
					// also need to look into those as well.
					for (int i=0;i<fast5files.length;i++) {
						if (fast5files[i].isDirectory()) {
							File [] subFast5files = fast5files[i].listFiles();
							
							for (int j=0;j<subFast5files.length;i++) {
								if (subFast5files[j].getName().endsWith(".fast5")) {
									files.add(subFast5files[j]);
								}
							}

						}
					}
					
				}
				else {
					files.add(file);
				}
			}
		}
		
		
		File [][] fileGroups;
		
		// See if we need to group together files from a casava group
		if (FastQCConfig.getInstance().casava) {
			fileGroups = CasavaBasename.getCasavaGroups(files.toArray(new File[0]));
		}
		else if (FastQCConfig.getInstance().nano) {
			fileGroups = NanoporeBasename.getNanoporeGroups(files.toArray(new File[0]));
		}
		else {
			fileGroups = new File [files.size()][1];
			for (int f=0;f<files.size();f++) {
				fileGroups[f][0] = files.elementAt(f);
			}
		}
		
		
		filesRemaining = new AtomicInteger(fileGroups.length);
				
		for (int i=0;i<fileGroups.length;i++) {

			try {
				processFile(fileGroups[i]);
			}
			catch (OutOfMemoryError e) {
				System.err.println("Ran out of memory for "+fileGroups[i][0]);
				e.printStackTrace();
				System.exit(2);
			}
			catch (Exception e) {
				System.err.println("Failed to process "+fileGroups[i][0]);
				e.printStackTrace();
				filesRemaining.decrementAndGet();
				somethingFailed = true;
			}
			

		}
		
		// We need to hold this class open as otherwise the main method
		// exits when it's finished.
		while (filesRemaining.intValue() > 0) {
			try {
				Thread.sleep(1000);
			} 
			catch (InterruptedException e) {}
		}
		if (somethingFailed) {
			System.exit(1);
		}
		System.exit(0);
		
	}
	
	public void processFile (File [] files) throws Exception {
		for (int f=0;f<files.length;f++) {
			if (!files[f].getName().startsWith("stdin") && !files[f].exists()) {
				throw new IOException(files[f].getName()+" doesn't exist");
			}
		}
		SequenceFile sequenceFile = SequenceFactory.getSequenceFile(files);			
						
		AnalysisRunner runner = new AnalysisRunner(sequenceFile);
		runner.addAnalysisListener(this);
			
		QCModule [] module_list = ModuleFactory.getStandardModuleList();

		runner.startAnalysis(module_list);

	}	
	
	public void analysisComplete(SequenceFile file, QCModule[] results) {
		File reportFile;
		
		if (showUpdates) System.out.println("Analysis complete for "+file.name());

		
		if (FastQCConfig.getInstance().output_dir != null) {
			String fileName = file.getFile().getName().replaceAll("stdin:","").replaceAll("\\.gz$","").replaceAll("\\.bz2$","").replaceAll("\\.txt$","").replaceAll("\\.fastq$", "").replaceAll("\\.fq$", "").replaceAll("\\.csfastq$", "").replaceAll("\\.sam$", "").replaceAll("\\.bam$", "")+"_fastqc.html";
			reportFile = new File(FastQCConfig.getInstance().output_dir+"/"+fileName);						
		}
		else {
			reportFile = new File(file.getFile().getAbsolutePath().replaceAll("stdin:","").replaceAll("\\.gz$","").replaceAll("\\.bz2$","").replaceAll("\\.txt$","").replaceAll("\\.fastq$", "").replaceAll("\\.fq$", "").replaceAll("\\.csfastq$", "").replaceAll("\\.sam$", "").replaceAll("\\.bam$", "")+"_fastqc.html");			
		}
		
		try {
			new HTMLReportArchive(file, results, reportFile);
		}
		catch (Exception e) {
			analysisExceptionReceived(file, e);
			return;
		}
		filesRemaining.decrementAndGet();

	}

	public void analysisUpdated(SequenceFile file, int sequencesProcessed, int percentComplete) {
		
		if (percentComplete % 5 == 0) {
			if (percentComplete == 105) {
				if (showUpdates) System.err.println("It seems our guess for the total number of records wasn't very good.  Sorry about that.");
			}
			if (percentComplete > 100) {
				if (showUpdates) System.err.println("Still going at "+percentComplete+"% complete for "+file.name());
			}
			else {
				if (showUpdates) System.err.println("Approx "+percentComplete+"% complete for "+file.name());
			}
		}
	}

	public void analysisExceptionReceived(SequenceFile file, Exception e) {
		System.err.println("Failed to process file "+file.name());
		somethingFailed = true;
		e.printStackTrace();
		filesRemaining.decrementAndGet();
	}

	public void analysisStarted(SequenceFile file) {
		if (showUpdates) System.err.println("Started analysis of "+file.name());
		
	}
	
}
