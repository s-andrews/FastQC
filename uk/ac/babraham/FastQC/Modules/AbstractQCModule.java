/**
 * Copyright Copyright 2010-12 Pierre Lindenbaum
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

package uk.ac.babraham.FastQC.Modules;

import java.awt.Graphics;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.zip.ZipEntry;
import java.util.zip.ZipOutputStream;

import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.imageio.ImageIO;
import javax.swing.JPanel;
import javax.swing.table.TableModel;

import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Report.HTMLReportArchive;
import uk.ac.babraham.FastQC.Utilities.ImageToBase64;
import uk.ac.babraham.FastQC.Utilities.ImageSaver.SVGImageSaver;

public abstract class AbstractQCModule implements QCModule {

	protected 	void simpleXhtmlReport(HTMLReportArchive report,String svgData, BufferedImage image, String alt) throws XMLStreamException {
		XMLStreamWriter xhtml = report.xhtmlStream();
		xhtml.writeStartElement("p");
		xhtml.writeEmptyElement("img");
		xhtml.writeAttribute("class", "indented");
		if (FastQCConfig.getInstance().svg_output) {
			xhtml.writeAttribute("src", ImageToBase64.svgImageToBase64(svgData));
		}
		else {
			xhtml.writeAttribute("src", ImageToBase64.imageToBase64(image));
		}
		xhtml.writeAttribute("alt", alt);
		
//		if(svgData!=null){
//			xhtml.writeAttribute("width",String.valueOf(img.getWidth()));
//			xhtml.writeAttribute("height",String.valueOf(img.getHeight()));
//		}
		
		xhtml.writeEndElement();//p
	}
	
	protected void writeDefaultImage (HTMLReportArchive report, String fileName, String imageTitle, int width, int height) throws IOException, XMLStreamException {
		ZipOutputStream zip = report.zipFile();
		
		// Write out the svg version of the image		
		JPanel resultsPanel = getResultsPanel();
		resultsPanel.setSize(width,height);
		resultsPanel.validate();

		
		String svgFilename = fileName.replaceAll("\\.png$", ".svg");
		zip.putNextEntry(new ZipEntry(report.folderName()+"/Images/"+svgFilename));
				
		String svgData = SVGImageSaver.saveImage(resultsPanel, zip);
		zip.closeEntry();

		// Write out the png version of the image
		zip.putNextEntry(new ZipEntry(report.folderName()+"/Images/"+fileName));
		BufferedImage b = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics g = b.createGraphics();

		resultsPanel.setDoubleBuffered(false);
		resultsPanel.addNotify();
		resultsPanel.validate();
		
		resultsPanel.print(g);

		g.dispose();
		
		ImageIO.write(b, "PNG", zip);
		zip.closeEntry();

		
		simpleXhtmlReport(report, svgData, b, imageTitle);

	}

	protected void writeSpecificImage (HTMLReportArchive report, JPanel resultsPanel, String fileName, String imageTitle, int width, int height) throws IOException, XMLStreamException {
		ZipOutputStream zip = report.zipFile();
		
		// Write out the svg version of the image		
		resultsPanel.setSize(width,height);
		resultsPanel.validate();

		
		String svgFilename = fileName.replaceAll("\\.png$", ".svg");
		zip.putNextEntry(new ZipEntry(report.folderName()+"/Images/"+svgFilename));
				
		String svgData = SVGImageSaver.saveImage(resultsPanel, zip);
		zip.closeEntry();

		// Write out the png version of the image
		zip.putNextEntry(new ZipEntry(report.folderName()+"/Images/"+fileName));
		BufferedImage b = new BufferedImage(width, height, BufferedImage.TYPE_INT_RGB);
		Graphics g = b.createGraphics();

		resultsPanel.setDoubleBuffered(false);
		resultsPanel.addNotify();
		resultsPanel.validate();
		
		resultsPanel.print(g);

		g.dispose();
		
		ImageIO.write(b, "PNG", zip);
		zip.closeEntry();

		
		simpleXhtmlReport(report, svgData, b, imageTitle);
	}


	protected void writeTable(HTMLReportArchive report, TableModel table) throws IOException,XMLStreamException {
		writeXhtmlTable(report,table);
		writeTextTable(report,table);	
	}

	protected void writeXhtmlTable(HTMLReportArchive report, TableModel table) throws IOException,XMLStreamException {
		XMLStreamWriter w=report.xhtmlStream();
		w.writeStartElement("table");
		w.writeStartElement("thead");
		w.writeStartElement("tr");
		
		for (int c=0;c<table.getColumnCount();c++) {
			w.writeStartElement("th");
			w.writeCharacters(table.getColumnName(c));
			w.writeEndElement();
		}
		
		w.writeEndElement();//tr
		w.writeEndElement();//thead
		w.writeStartElement("tbody");
		
		for (int r=0;r<table.getRowCount();r++) {
			w.writeStartElement("tr");
			for (int c=0;c<table.getColumnCount();c++) {
				w.writeStartElement("td");
				w.writeCharacters(String.valueOf(table.getValueAt(r, c)));
				w.writeEndElement();//td
			}
			w.writeEndElement();//tr
		}
		w.writeEndElement();//tbody
		w.writeEndElement();
	}


	protected void writeTextTable(HTMLReportArchive report, TableModel table) throws IOException {
		StringBuffer d=report.dataDocument();
		d.append("#");

		for (int c=0;c<table.getColumnCount();c++) {
			if (c!=0) d.append("\t");
			d.append(table.getColumnName(c));
		}
		
		d.append("\n");

		// Do the rows
		for (int r=0;r<table.getRowCount();r++) {
			for (int c=0;c<table.getColumnCount();c++) {
				if (c!=0) d.append("\t");
				d.append(table.getValueAt(r, c));
			}
			d.append("\n");
		}

	}

}
