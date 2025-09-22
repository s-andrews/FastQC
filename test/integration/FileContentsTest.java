package test.integration;

import org.approvaltests.Approvals;
import org.approvaltests.core.Options;
import org.approvaltests.core.Scrubber;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static test.integration.TestHelpers.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import test.integration.HtmlContentHelpers;
import test.integration.cli.Cli;
import test.integration.cli.CliScenario;
public class FileContentsTest {

    @ParameterizedTest
    @MethodSource("test.integration.cli.CliScenario#scenarios")
    public void verify_data(String name) throws Exception {
        var content = ExecuteAndExtractFileContent(name, "fastqc_data.txt");

        var options = new Options()
            .forFile()
            .withName("FileContentsTest_" + name + "_fastqc_data", " txt"); 
        Approvals.verify(content, options);
    }

    @ParameterizedTest
    @MethodSource("test.integration.cli.CliScenario#scenarios")
    public void verify_html(String name) throws Exception {
        var content = ExecuteAndExtractFileContent(name, "fastqc_report.html");
        content = HtmlContentHelpers.normalizeHtml(content);
        
        var options = new Options()
            .withScrubber(HtmlContentHelpers::scrubDates)
            .forFile()
            .withName("FileContentsTest_" + name + "_fastqc_report", " html"); 
        Approvals.verify(content, options);
    }

    private String ExecuteAndExtractFileContent(String scenarioName, String extractFileName)
            throws Exception, IOException {
        var outputDir = new File(CliScenario.TEST_OUT_DIR).toPath();
        if (!Files.exists(outputDir)) {
            Files.createDirectories(outputDir);
        }

        // create the scenario
        var scenario = new CliScenario(scenarioName, new String[] {
                "fastqc.output_dir=" + outputDir
        });

        // execute the scenario, checking its output and that it is successful
        Cli.Execute(scenario).assertSuccess();
        TestHelpers.unzip(scenario.ZipFilePath, scenario.ZipExtractionPath);

        // locate the inner folder and check it exists
        var folder = new File(scenario.ZipExtractionPath, scenario.ZipInnerFolderName).getAbsolutePath().toString();
        assertFolderExists(folder);

        // locate the data file
        var extractedFile = new File(folder, extractFileName);
        assertFileExists(extractedFile.getAbsolutePath());

        System.out.println("Extracted file: " + extractedFile);
        return Files.readString(extractedFile.toPath());
    }
}