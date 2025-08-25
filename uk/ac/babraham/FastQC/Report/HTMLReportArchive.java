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
import java.util.HashMap;
import java.util.Map;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;
import java.util.zip.ZipOutputStream;

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

public class HTMLReportArchive {
	private XMLStreamWriter xhtml=null;
	private StringBuffer data = new StringBuffer();
	private QCModule [] modules;
	private ZipOutputStream zip;
	private SequenceFile sequenceFile;
	private byte [] buffer = new byte[1024];
	private File htmlFile;
	private File zipFile;
	private StringWriter moduleContentWriter;
	private String htmlTemplate;
	private Map<String, String> helpFileMapping;

	public HTMLReportArchive (SequenceFile sequenceFile, QCModule [] modules, File htmlFile) throws IOException, XMLStreamException {
		this.sequenceFile = sequenceFile;
		this.modules = modules;
		this.htmlFile = htmlFile;
		this.zipFile = new File(htmlFile.getAbsoluteFile().toString().replaceAll("\\.html$", "")+".zip");

		// Load the HTML template
		this.htmlTemplate = loadTemplate("/Templates/report_template.html");

		// Initialize help file mapping
		this.helpFileMapping = initializeHelpFileMapping();

		// Create XMLStreamWriter for module content
		this.moduleContentWriter = new StringWriter();
		XMLOutputFactory xmlfactory = XMLOutputFactory.newInstance();
		this.xhtml= xmlfactory.createXMLStreamWriter(moduleContentWriter);


		zip = new ZipOutputStream(new FileOutputStream(zipFile));
		zip.putNextEntry(new ZipEntry(folderName()+"/"));
		zip.putNextEntry(new ZipEntry(folderName()+"/Icons/"));
		zip.putNextEntry(new ZipEntry(folderName()+"/Images/"));

		// Initialize data document
		data.append("##FastQC\t");
		data.append(FastQCApplication.VERSION);
		data.append("\n");

		// Generate module content using template-based approach
		String moduleContent = generateModuleContent();

		// Close XMLStreamWriter (no longer needed for main content)
		xhtml.flush();
		xhtml.close();

		// Generate final HTML from template
		String finalHtml = generateHtmlFromTemplate(moduleContent);

		zip.putNextEntry(new ZipEntry(folderName()+"/fastqc_report.html"));
		zip.write(finalHtml.getBytes());
		zip.closeEntry();
		zip.putNextEntry(new ZipEntry(folderName()+"/fastqc_data.txt"));
		zip.write(data.toString().getBytes());
		zip.closeEntry();

		// Generate summary.txt
		String summaryText = generateSummaryText();
		zip.putNextEntry(new ZipEntry(folderName()+"/summary.txt"));
		zip.write(summaryText.getBytes());
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

	private String loadTemplate(String templatePath) throws IOException {
		InputStream templateStream = getClass().getResourceAsStream(templatePath);
		if (templateStream == null) {
			throw new IOException("Template not found: " + templatePath);
		}

		StringWriter templateWriter = new StringWriter();
		byte[] buffer = new byte[1024];
		int nRead;
		while ((nRead = templateStream.read(buffer)) != -1) {
			templateWriter.write(new String(buffer, 0, nRead));
		}
		templateStream.close();
		return templateWriter.toString();
	}

		private String loadCss() throws IOException {
		return loadTemplate("/Templates/fastqc.css");
	}

	private String generateSummaryItems() throws IOException {
		StringBuffer summaryItems = new StringBuffer();
		String sidebarItemTemplate = loadTemplate("/Templates/sidebar_item.html");

		for (int m=0;m<modules.length;m++) {
			if (modules[m].ignoreInReport()) {
				continue;
			}

			String item = sidebarItemTemplate;
			item = item.replace("{{MODULE_INDEX}}", String.valueOf(m));
			item = item.replace("{{MODULE_NAME}}", modules[m].name());

			if (modules[m].raisesError()) {
				item = item.replace("{{STATUS_CLASS}}", "sidebar-error");
				item = item.replace("{{STATUS_TEXT}}", "Error");
			} else if (modules[m].raisesWarning()) {
				item = item.replace("{{STATUS_CLASS}}", "sidebar-warning");
				item = item.replace("{{STATUS_TEXT}}", "Warn");
			} else {
				item = item.replace("{{STATUS_CLASS}}", "sidebar-pass");
				item = item.replace("{{STATUS_TEXT}}", "Pass");
			}

			summaryItems.append(item).append("\n");
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
			} else if (modules[m].raisesWarning()) {
				summaryText.append("WARN");
			} else {
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

		private String getStatusIcon(QCModule module) throws IOException {
		if (module.raisesError()) {
			return loadTemplate("/Templates/Icons/error.svg");
		} else if (module.raisesWarning()) {
			return loadTemplate("/Templates/Icons/warning.svg");
		} else {
			return loadTemplate("/Templates/Icons/pass.svg");
		}
	}

	private String getFastQCIconWithUniqueIds(String suffix) throws IOException {
		String svgContent = loadTemplate("/Templates/Icons/fastqc_icon.svg");

		// Replace all IDs and their references with unique versions
		return svgContent.replaceAll("id=\"([^\"]+)\"", "id=\"$1_" + suffix + "\"")
				         .replaceAll("url\\(#([^)]+)\\)", "url(#$1_" + suffix + ")");
	}

	private String generateModuleContent() throws IOException, XMLStreamException {
		StringBuffer allModulesContent = new StringBuffer();
		String moduleWrapperTemplate = loadTemplate("/Templates/module_wrapper.html");

		// Generate content for each module
		for (int m=0;m<modules.length;m++) {
			if (modules[m].ignoreInReport()) continue;

			// Create a separate XMLStreamWriter for this module's content
			StringWriter moduleBodyWriter = new StringWriter();
			XMLOutputFactory xmlfactory = XMLOutputFactory.newInstance();
			XMLStreamWriter moduleXhtml = xmlfactory.createXMLStreamWriter(moduleBodyWriter);

			// Temporarily switch the xhtml writer for this module
			XMLStreamWriter originalXhtml = this.xhtml;
			this.xhtml = moduleXhtml;

			// Add data for this module
			data.append(">>");
			data.append(modules[m].name());
			data.append("\t");
			if (modules[m].raisesError()) {
				data.append("fail");
			} else if (modules[m].raisesWarning()) {
				data.append("warn");
			} else {
				data.append("pass");
			}
			data.append("\n");

			// Let the module generate its content
			modules[m].makeReport(this);
			data.append(">>END_MODULE\n");

			// Close the module's XMLStreamWriter
			moduleXhtml.flush();
			moduleXhtml.close();

			// Restore the original XMLStreamWriter
			this.xhtml = originalXhtml;

			// Apply the module wrapper template
			String moduleWrapper = moduleWrapperTemplate;
			moduleWrapper = moduleWrapper.replace("{{MODULE_INDEX}}", String.valueOf(m));
			moduleWrapper = moduleWrapper.replace("{{MODULE_NAME}}", modules[m].name());
			moduleWrapper = moduleWrapper.replace("{{STATUS_ICON}}", getStatusIcon(modules[m]));
			moduleWrapper = moduleWrapper.replace("{{HELP_CONTENT}}", extractHelpText(modules[m].name()));
			moduleWrapper = moduleWrapper.replace("{{MODULE_CONTENT}}", moduleBodyWriter.toString());

			allModulesContent.append(moduleWrapper).append("\n");
		}

		return allModulesContent.toString();
	}

	private String generateHtmlFromTemplate(String moduleContent) throws IOException {
		SimpleDateFormat df = new SimpleDateFormat("EEE d MMM yyyy");

		String html = htmlTemplate;
		html = html.replace("{{TITLE}}", sequenceFile.name() + " FastQC Report");
		html = html.replace("{{CSS_CONTENT}}", loadCss());
		html = html.replace("{{DATE}}", df.format(new Date()));
		html = html.replace("{{FILENAME}}", sequenceFile.name());
		html = html.replace("{{FASTQC_ICON_SVG_MOBILE}}", getFastQCIconWithUniqueIds("mobile"));
		html = html.replace("{{FASTQC_ICON_SVG_SIDEBAR}}", getFastQCIconWithUniqueIds("sidebar"));
		html = html.replace("{{SUMMARY_ITEMS}}", generateSummaryItems());
		html = html.replace("{{MODULE_CONTENT}}", moduleContent);
		html = html.replace("{{VERSION}}", FastQCApplication.VERSION);

		return html;
	}






	/**
	 * Initialize mapping between module names and their corresponding help file paths
	 */
	private Map<String, String> initializeHelpFileMapping() {
		Map<String, String> mapping = new HashMap<String, String>();
		mapping.put("Basic statistics", "/Help/3 Analysis Modules/1 Basic statistics.html");
		mapping.put("Per base sequence quality", "/Help/3 Analysis Modules/2 Per Base Sequence Quality.html");
		mapping.put("Per sequence quality scores", "/Help/3 Analysis Modules/3 Per Sequence Quality Scores.html");
		mapping.put("Per base sequence content", "/Help/3 Analysis Modules/4 Per Base Sequence Content.html");
		mapping.put("Per sequence GC content", "/Help/3 Analysis Modules/5 Per Sequence GC Content.html");
		mapping.put("Per base N content", "/Help/3 Analysis Modules/6 Per Base N Content.html");
		mapping.put("Sequence length distribution", "/Help/3 Analysis Modules/7 Sequence length distribution.html");
		mapping.put("Sequence duplication levels", "/Help/3 Analysis Modules/8 Duplicate Sequences.html");
		mapping.put("Overrepresented sequences", "/Help/3 Analysis Modules/9 Overrepresented Sequences.html");
		mapping.put("Adapter content", "/Help/3 Analysis Modules/10 Adapter content.html");
		mapping.put("Kmer Content", "/Help/3 Analysis Modules/11 Kmer Content.html");
		mapping.put("Per tile sequence quality", "/Help/3 Analysis Modules/12 Per Tile Sequence Quality.html");
		return mapping;
	}

	/**
	 * Extract text content from help HTML files, removing images and HTML tags
	 */
	private String extractHelpText(String moduleName) {
		String helpFilePath = helpFileMapping.get(moduleName);
		if (helpFilePath == null) {
			return "<p>Help documentation not available for this module.</p>";
		}

		try {
			String helpHtml = loadTemplate(helpFilePath);
			return convertHelpHtmlToText(helpHtml);
		} catch (IOException e) {
			return "<p>Help documentation could not be loaded for this module.</p>";
		}
	}

	/**
	 * Convert help HTML to clean text content, removing images and preserving structure
	 */
	private String convertHelpHtmlToText(String htmlContent) {
		// Remove the HTML document structure, head, and body tags
		String content = htmlContent.replaceAll("(?s)<html.*?>", "");
		content = content.replaceAll("(?s)<head.*?>.*?</head>", "");
		content = content.replaceAll("(?s)<body[^>]*>", "");
		content = content.replaceAll("(?s)</body>.*?</html>", "");

		// Remove all img tags and p tags containing only images
		content = content.replaceAll("(?s)<p>\\s*<img[^>]*>\\s*</p>", "");
		content = content.replaceAll("(?s)<img[^>]*>", "");

		// Remove empty p tags
		content = content.replaceAll("(?s)<p>\\s*</p>", "");

		// Remove h1 tags completely since they duplicate the module name
		content = content.replaceAll("(?s)<h1[^>]*>.*?</h1>", "");

		// Convert h2 tags to h4 for better integration
		content = content.replaceAll("(?s)<h2([^>]*)>", "<h4$1>");
		content = content.replaceAll("(?s)</h2>", "</h4>");

		// Clean up extra whitespace and newlines
		content = content.replaceAll("\\n\\s*\\n\\s*\\n", "\n\n");
		content = content.replaceAll("(?m)^\\s+", "");
		content = content.trim();

		return content;
	}

	/**
	 * Ensure proper HTML list structure by closing any unclosed li tags
	 */
	private String ensureProperListStructure(String content) {
		// Split content into lines for processing
		String[] lines = content.split("\n");
		StringBuilder result = new StringBuilder();
		boolean inListItem = false;

		for (String line : lines) {
			String trimmed = line.trim();

			// Check if this line starts a new list item
			if (trimmed.startsWith("<li>") || trimmed.matches("^<li[^>]*>.*")) {
				// Close previous list item if open
				if (inListItem) {
					result.append("</li>\n");
				}
				result.append(line).append("\n");
				inListItem = true;
			}
			// Check if this line closes a list item
			else if (trimmed.equals("</li>")) {
				result.append(line).append("\n");
				inListItem = false;
			}
			// Check if this line ends a list
			else if (trimmed.equals("</ol>") || trimmed.equals("</ul>")) {
				// Close current list item if open
				if (inListItem) {
					result.append("</li>\n");
					inListItem = false;
				}
				result.append(line).append("\n");
			}
			// Regular line
			else {
				result.append(line).append("\n");
			}
		}

		// Close final list item if still open
		if (inListItem) {
			result.append("</li>\n");
		}

		return result.toString();
	}

}
