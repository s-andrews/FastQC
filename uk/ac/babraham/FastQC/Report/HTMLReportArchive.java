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
package uk.ac.babraham.FastQC.Report;

import java.awt.image.BufferedImage;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.PrintWriter;
import java.io.StringReader;
import java.io.StringWriter;
import java.text.SimpleDateFormat;
import java.util.Date;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

import javax.imageio.ImageIO;
import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.stream.XMLOutputFactory;
import javax.xml.stream.XMLStreamException;
import javax.xml.stream.XMLStreamWriter;
import javax.xml.transform.Templates;
import javax.xml.transform.TransformerFactory;
import javax.xml.transform.dom.DOMSource;
import javax.xml.transform.stream.StreamResult;

import org.w3c.dom.Document;
import org.xml.sax.InputSource;

import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.FastQCApplication;
import uk.ac.babraham.FastQC.Modules.QCModule;
import uk.ac.babraham.FastQC.Sequence.SequenceFile;
import uk.ac.babraham.FastQC.Utilities.ImageToBase64;

public class HTMLReportArchive {
	private XMLStreamWriter xhtml=null;
	private StringBuffer data = new StringBuffer();
	private QCModule [] modules;
	private ZipOutputStream zip;
	private SequenceFile sequenceFile;
	private byte [] buffer = new byte[1024];
	private File htmlFile;
	private File zipFile;

	public HTMLReportArchive (SequenceFile sequenceFile, QCModule [] modules, File htmlFile) throws IOException, XMLStreamException {
		this.sequenceFile = sequenceFile;
		this.modules = modules;
		this.htmlFile = htmlFile;
		this.zipFile = new File(htmlFile.getAbsoluteFile().toString().replaceAll("\\.html$", "")+".zip");
		StringWriter htmlStr = new StringWriter();
		XMLOutputFactory xmlfactory = XMLOutputFactory.newInstance();
		this.xhtml= xmlfactory.createXMLStreamWriter(htmlStr);
		
		
		zip = new ZipOutputStream(new FileOutputStream(zipFile));
		zip.putNextEntry(new ZipEntry(folderName()+"/"));
		zip.putNextEntry(new ZipEntry(folderName()+"/Icons/"));
		zip.putNextEntry(new ZipEntry(folderName()+"/Images/"));
		startDocument();
		for (int m=0;m<modules.length;m++) {
			
			if (modules[m].ignoreInReport()) continue;
			
			xhtml.writeStartElement("div");
			xhtml.writeAttribute("class", "module");
			xhtml.writeStartElement("h2");
			xhtml.writeAttribute("id", "M"+m);
			
			
			// Add an icon before the module name
			xhtml.writeStartElement("svg");
			xhtml.writeAttribute("xmlns", "http://www.w3.org/2000/svg");
			xhtml.writeAttribute("width", "24");
			xhtml.writeAttribute("height", "24");
			xhtml.writeAttribute("viewBox", "0 0 24 24");
			xhtml.writeStartElement("path");
			if (modules[m].raisesError())
				{
				xhtml.writeAttribute("fill", "#d65d3e");
				xhtml.writeAttribute("d", "M12 2c5.53 0 10 4.47 10 10s-4.47 10-10 10S2 17.53 2 12S6.47 2 12 2m3.59 5L12 10.59L8.41 7L7 8.41L10.59 12L7 15.59L8.41 17L12 13.41L15.59 17L17 15.59L13.41 12L17 8.41z");
				}

			else if (modules[m].raisesWarning())
				{
				xhtml.writeAttribute("fill", "#eab30d");
				xhtml.writeAttribute("d", "M13 14h-2V9h2m0 9h-2v-2h2M1 21h22L12 2z");
				}
			else {
				xhtml.writeAttribute("fill", "#4ba359");
				xhtml.writeAttribute("d", "M12 2C6.5 2 2 6.5 2 12s4.5 10 10 10s10-4.5 10-10S17.5 2 12 2m-2 15l-5-5l1.41-1.41L10 14.17l7.59-7.59L19 8z");
				}
			xhtml.writeEndElement(); // path
			xhtml.writeEndElement(); // svg

			
			xhtml.writeCharacters(modules[m].name());
			data.append(">>");
			data.append(modules[m].name());
			data.append("\t");
			if (modules[m].raisesError()) {
				data.append("fail");
			}
			else if (modules[m].raisesWarning()) {
				data.append("warn");
			}
			else {
				data.append("pass");
			}
			data.append("\n");
			xhtml.writeEndElement();
			modules[m].makeReport(this);
			data.append(">>END_MODULE\n");
			
			xhtml.writeEndElement();
		}
		closeDocument();
		
		zip.putNextEntry(new ZipEntry(folderName()+"/fastqc_report.html"));
		xhtml.flush();
		xhtml.close();
		zip.write(htmlStr.toString().getBytes());
		zip.closeEntry();
		zip.putNextEntry(new ZipEntry(folderName()+"/fastqc_data.txt"));
		zip.write(data.toString().getBytes());
		zip.closeEntry();
		
		//XSL-FO
		try {
			DocumentBuilderFactory domFactory=DocumentBuilderFactory.newInstance();
			domFactory.setNamespaceAware(false);
			DocumentBuilder builder=domFactory.newDocumentBuilder();
			Document src=builder.parse(new InputSource( new StringReader(htmlStr.toString())));
			InputStream rsrc=getClass().getResourceAsStream("/Templates/fastqc2fo.xsl");
			if(rsrc!=null)
				{
				domFactory.setNamespaceAware(true);
				builder=domFactory.newDocumentBuilder();
				Document html2fo=builder.parse(rsrc);
				rsrc.close();
				
				TransformerFactory tf=TransformerFactory.newInstance();
				Templates templates=tf.newTemplates(new DOMSource(html2fo));
				zip.putNextEntry(new ZipEntry(folderName()+"/fastqc.fo"));
				templates.newTransformer().transform(new DOMSource(src), new StreamResult(zip));
				zip.closeEntry();
				}
			}
		catch (Exception e) {
			e.printStackTrace();
		}
		
		
		zip.close();
		
		// Save the HTML file at the same level as the zip file
		
		PrintWriter pr = new PrintWriter(new FileWriter(htmlFile));
		
		pr.print(htmlStr.toString());
		
		pr.close();

		
		if (FastQCConfig.getInstance().do_unzip) {
			unzipZipFile(zipFile);
			if (FastQCConfig.getInstance().delete_after_unzip) {
				zipFile.delete();
			}
		}
	}
	
	private void unzipZipFile (File file) throws IOException {
		ZipFile zipFile = new ZipFile(file);
		Enumeration<? extends ZipEntry> entries = zipFile.entries();
		int size;
		byte [] buffer = new byte[1024];
		
		while (entries.hasMoreElements()) {
			ZipEntry entry = entries.nextElement();
			
//			System.out.println("Going to extract '"+entry.getName()+"'");
			
			if (entry.isDirectory()) {
				File dir = new File(file.getParent()+"/"+entry.getName());
				if (dir.exists() && dir.isDirectory()) continue; // Don't need to do anything
				if (dir.exists() && ! dir.isDirectory()) throw new IOException ("File exists with dir name "+dir.getName());
				if (!dir.mkdir()) throw new IOException("Failed to make dir for "+dir.getName());
				continue;
			}

			BufferedInputStream bis = new BufferedInputStream(zipFile.getInputStream(entry));
			BufferedOutputStream bos = new BufferedOutputStream(new FileOutputStream(file.getParent()+"/"+entry.getName()),buffer.length);
			while ((size = bis.read(buffer,0,buffer.length)) != -1) {
				bos.write(buffer,0,size);
			}
			bos.flush();
			bos.close();
			bis.close();
		}
		
		zipFile.close();
	}

	public XMLStreamWriter xhtmlStream ()
		{
		return this.xhtml;
		}
	
	public StringBuffer dataDocument() {
		return data;
	}
	
	public String folderName () {
		return htmlFile.getName().replaceAll("\\.html$", "");
	}
	
	public ZipOutputStream zipFile () {
		return zip;
	}
	
	private void startDocument () throws IOException,XMLStreamException
		{
		
		// Just put the fastQC version at the start of the text report
		data.append("##FastQC\t");
		data.append(FastQCApplication.VERSION);
		data.append("\n");
		
		// Add in the icon files for pass/fail/warn
		for(String icnName:new String[]{
				"fastqc_icon.png"})
			{
			InputStream in =getClass().getResourceAsStream("/Templates/Icons/"+icnName);
			if(in==null) continue;
			zip.putNextEntry(new ZipEntry(folderName()+"/Icons/"+icnName));
			int len;
			while ((len = in.read(buffer)) > 0) {
				zip.write(buffer, 0, len);
			}
			in.close();
			zip.closeEntry();
			}
				
		

		SimpleDateFormat df = new SimpleDateFormat("EEE d MMM yyyy");
		
		xhtml.writeDTD("<!DOCTYPE html>");
		xhtml.writeStartElement("html");
		xhtml.writeStartElement("head");
		
		xhtml.writeStartElement("title");
		xhtml.writeCharacters(sequenceFile.name());
		xhtml.writeCharacters(" FastQC Report");
		xhtml.writeEndElement();//title
		
		InputStream rsrc=getClass().getResourceAsStream("/Templates/header_template.html");
		if(rsrc!=null)
			{
			xhtml.writeStartElement("style");
			xhtml.writeAttribute("type", "text/css");

			byte array[]=new byte[128];
			int nRead;
			while((nRead=rsrc.read(array))!=-1) { xhtml.writeCharacters(new String(array,0,nRead));}
			rsrc.close();
			xhtml.writeEndElement();//style
			}		

		
		
		
		xhtml.writeEndElement();//head
		
		xhtml.writeStartElement("body");
		
		xhtml.wrxiteStartElement("div");
		xhtml.writeAttribute("class", "header");
		
		xhtml.writeStartElement("div");
		xhtml.writeAttribute("id", "header_title");
		
		xhtml.writeEmptyElement("img");
		xhtml.writeAttribute("src", base64ForIcon("Icons/fastqc_icon.png"));
		xhtml.writeAttribute("alt", "FastQC");
		xhtml.writeCharacters("FastQC Report");
		xhtml.writeEndElement();//div
		
		xhtml.writeStartElement("div");
		xhtml.writeAttribute("id", "header_filename");
		xhtml.writeStartElement("p");
		xhtml.writeCharacters(df.format(new Date()));
		xhtml.writeEndElement();//p
		xhtml.writeStartElement("p");
		xhtml.writeStartElement("strong");
		xhtml.writeCharacters(sequenceFile.name());
		xhtml.writeEndElement();//strong
		xhtml.writeEndElement();//p
		xhtml.writeEndElement();//div
		xhtml.writeEndElement();//div
		
		
		xhtml.writeStartElement("div");
		xhtml.writeAttribute("class", "summary");
		
		xhtml.writeStartElement("h2");
		xhtml.writeCharacters("Summary");
		xhtml.writeEndElement();//h2
		
		
		xhtml.writeStartElement("ul");
		
		StringBuffer summaryText = new StringBuffer();
		
		for (int m=0;m<modules.length;m++) {
			
			if (modules[m].ignoreInReport()) {
				continue;
			}
			xhtml.writeStartElement("li");
			xhtml.writeStartElement("a");
			xhtml.writeAttribute("href", "#M"+m);

			xhtml.writeStartElement("span");
			if (modules[m].raisesError()) {
				xhtml.writeAttribute("class", "sidebar-error");
				xhtml.writeCharacters("Error");
				summaryText.append("FAIL");
				}
			else if (modules[m].raisesWarning()) {
				xhtml.writeAttribute("class", "sidebar-warning");
				xhtml.writeCharacters("Warn");
				summaryText.append("WARN");
			}
			else {
				xhtml.writeAttribute("class", "sidebar-pass");
				xhtml.writeCharacters("Pass");
				summaryText.append("PASS");
			}
			xhtml.writeEndElement(); // span

			xhtml.writeCharacters(modules[m].name());
			xhtml.writeEndElement();//a
			xhtml.writeEndElement();//li

			summaryText.append("\t");
			summaryText.append(modules[m].name());
			summaryText.append("\t");
			summaryText.append(sequenceFile.name());
			summaryText.append(FastQCConfig.getInstance().lineSeparator);


		}
		xhtml.writeEndElement();//ul
		xhtml.writeEndElement();//div
		
		xhtml.writeStartElement("div");
		xhtml.writeAttribute("class", "main");
		

		zip.putNextEntry(new ZipEntry(folderName()+"/summary.txt"));
		zip.write(summaryText.toString().getBytes());

	}
	
	private String base64ForIcon (String path) {
		try {
			BufferedImage b = ImageIO.read(ClassLoader.getSystemResource("Templates/"+path));
			return (ImageToBase64.imageToBase64(b));
		}
		catch (IOException ioe) {
			ioe.printStackTrace();
			return "Failed";
		}
	}
	
	private void closeDocument () throws XMLStreamException
		{
		xhtml.writeEndElement();//div
		xhtml.writeStartElement("div");
		xhtml.writeAttribute("class", "footer");
		xhtml.writeCharacters("Produced by ");
		xhtml.writeStartElement("a");
		xhtml.writeAttribute("href", "http://www.bioinformatics.babraham.ac.uk/projects/fastqc/");
		xhtml.writeCharacters("FastQC");
		xhtml.writeEndElement();//a
		xhtml.writeCharacters("  (version "+FastQCApplication.VERSION+")");
		xhtml.writeEndElement();//div
		
		xhtml.writeEndElement();//body
		xhtml.writeEndElement();//html
		}
	
	
	
	
}