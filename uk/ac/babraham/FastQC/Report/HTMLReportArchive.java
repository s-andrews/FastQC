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
	private final String reportTemplate;
	private final String sidebarItemTemplate;
	private final String moduleWrapperTemplate;
	private final String cssContent;

	public HTMLReportArchive (SequenceFile sequenceFile, QCModule [] modules, File htmlFile) throws IOException, XMLStreamException {
		this.sequenceFile = sequenceFile;
		this.modules = modules;
		this.htmlFile = htmlFile;
		this.zipFile = new File(htmlFile.getAbsoluteFile().toString().replaceAll("\\.html$", "")+".zip");
		this.reportTemplate = loadTemplate("/Templates/report_template.html");
		this.sidebarItemTemplate = loadTemplate("/Templates/sidebar_item.html");
		this.moduleWrapperTemplate = loadTemplate("/Templates/module_wrapper.html");
		this.cssContent = loadResource("/Templates/fastqc.css");

		XMLOutputFactory xmlfactory = XMLOutputFactory.newInstance();
		this.xhtml= xmlfactory.createXMLStreamWriter(new StringWriter());
		
		
		zip = new ZipOutputStream(new FileOutputStream(zipFile));
		zip.putNextEntry(new ZipEntry(folderName()+"/"));
		zip.putNextEntry(new ZipEntry(folderName()+"/Icons/"));
		zip.putNextEntry(new ZipEntry(folderName()+"/Images/"));

		data.append("##FastQC\t");
		data.append(FastQCApplication.VERSION);
		data.append("\n");
		addIconsToZip();

		String moduleContent = generateModuleContent();
		xhtml.flush();
		xhtml.close();

		String finalHtml = generateHtmlFromTemplate(moduleContent);

		zip.putNextEntry(new ZipEntry(folderName()+"/fastqc_report.html"));
		zip.write(finalHtml.getBytes());
		zip.closeEntry();
		zip.putNextEntry(new ZipEntry(folderName()+"/fastqc_data.txt"));
		zip.write(data.toString().getBytes());
		zip.closeEntry();

		zip.putNextEntry(new ZipEntry(folderName()+"/summary.txt"));
		zip.write(generateSummaryText().getBytes());
		zip.closeEntry();

		//XSL-FO
		try {
			DocumentBuilderFactory domFactory=DocumentBuilderFactory.newInstance();
			domFactory.setNamespaceAware(false);
			DocumentBuilder builder=domFactory.newDocumentBuilder();
			Document src=builder.parse(new InputSource( new StringReader(finalHtml)));
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
		
		pr.print(finalHtml);
		
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
	
	private void addIconsToZip() throws IOException {
		for(String icnName:new String[]{
				"fastqc_icon.png",
				"warning.png",
				"error.png",
				"tick.png"})
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
	}

	private String loadResource(String resourcePath) throws IOException {
		InputStream in = getClass().getResourceAsStream(resourcePath);
		if (in == null) {
			throw new IOException("Resource not found: " + resourcePath);
		}

		StringWriter writer = new StringWriter();
		byte[] readBuffer = new byte[1024];
		int nRead;
		while ((nRead = in.read(readBuffer)) != -1) {
			writer.write(new String(readBuffer, 0, nRead));
		}
		in.close();
		return writer.toString();
	}

	private String loadTemplate(String templatePath) throws IOException {
		return loadResource(templatePath).stripTrailing();
	}

	private String escapeXml(String s) {
		return s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;");
	}

	private String generateSummaryItems() throws IOException {
		StringBuffer summaryItems = new StringBuffer();

		for (int m=0;m<modules.length;m++) {
			if (modules[m].ignoreInReport()) {
				continue;
			}

			String icon;
			String alt;
			if (modules[m].raisesError()) {
				icon = base64ForIcon("Icons/error.png");
				alt = "[FAIL]";
			}
			else if (modules[m].raisesWarning()) {
				icon = base64ForIcon("Icons/warning.png");
				alt = "[WARNING]";
			}
			else {
				icon = base64ForIcon("Icons/tick.png");
				alt = "[PASS]";
			}

			String item = sidebarItemTemplate
				.replace("{{MODULE_INDEX}}", String.valueOf(m))
				.replace("{{ICON}}", icon)
				.replace("{{ALT}}", alt)
				.replace("{{MODULE_NAME}}", escapeXml(modules[m].name()));

			summaryItems.append(item);
		}

		return summaryItems.toString();
	}

	private String generateSummaryText() {
		StringBuffer summaryText = new StringBuffer();

		for (int m=0;m<modules.length;m++) {
			if (modules[m].ignoreInReport()) {
				continue;
			}

			if (modules[m].raisesError()) {
				summaryText.append("FAIL");
			}
			else if (modules[m].raisesWarning()) {
				summaryText.append("WARN");
			}
			else {
				summaryText.append("PASS");
			}

			summaryText.append("\t");
			summaryText.append(modules[m].name());
			summaryText.append("\t");
			summaryText.append(sequenceFile.name());
			summaryText.append(FastQCConfig.getInstance().lineSeparator);
		}

		return summaryText.toString();
	}

	private String generateModuleContent() throws IOException, XMLStreamException {
		StringBuffer allModulesContent = new StringBuffer();

		for (int m=0;m<modules.length;m++) {
			if (modules[m].ignoreInReport()) continue;

			StringWriter moduleBodyWriter = new StringWriter();
			XMLOutputFactory xmlfactory = XMLOutputFactory.newInstance();
			XMLStreamWriter moduleXhtml = xmlfactory.createXMLStreamWriter(moduleBodyWriter);

			XMLStreamWriter originalXhtml = this.xhtml;
			this.xhtml = moduleXhtml;

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

			modules[m].makeReport(this);
			data.append(">>END_MODULE\n");

			moduleXhtml.flush();
			moduleXhtml.close();

			this.xhtml = originalXhtml;

			String icon;
			String alt;
			if (modules[m].raisesError()) {
				icon = base64ForIcon("Icons/error.png");
				alt = "[FAIL]";
			}
			else if (modules[m].raisesWarning()) {
				icon = base64ForIcon("Icons/warning.png");
				alt = "[WARN]";
			}
			else {
				icon = base64ForIcon("Icons/tick.png");
				alt = "[OK]";
			}

			String moduleWrapper = moduleWrapperTemplate
				.replace("{{MODULE_INDEX}}", String.valueOf(m))
				.replace("{{ICON}}", icon)
				.replace("{{ALT}}", alt)
				.replace("{{MODULE_NAME}}", escapeXml(modules[m].name()))
				.replace("{{MODULE_CONTENT}}", moduleBodyWriter.toString());

			allModulesContent.append(moduleWrapper);
		}

		return allModulesContent.toString();
	}

	private String generateHtmlFromTemplate(String moduleContent) throws IOException {
		SimpleDateFormat df = new SimpleDateFormat("EEE d MMM yyyy");

		return reportTemplate
			.replace("{{TITLE}}", escapeXml(sequenceFile.name()))
			.replace("{{CSS_CONTENT}}", escapeXml(cssContent))
			.replace("{{FASTQC_ICON}}", base64ForIcon("Icons/fastqc_icon.png"))
			.replace("{{DATE}}", escapeXml(df.format(new Date())))
			.replace("{{FILENAME}}", escapeXml(sequenceFile.name()))
			.replace("{{SUMMARY_ITEMS}}", generateSummaryItems())
			.replace("{{MODULE_CONTENT}}", moduleContent)
			.replace("{{VERSION}}", FastQCApplication.VERSION);
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
}
