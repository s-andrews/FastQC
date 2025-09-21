package test.integration;

import static org.junit.jupiter.api.Assertions.assertTrue;
import static test.integration.TestHelpers.*;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.stream.Stream;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import test.integration.cli.Cli;
import test.integration.cli.CliResult;
import test.integration.cli.CliScenario;

import static org.junit.jupiter.api.Assertions.*;

public class FolderStructureTest {

    @ParameterizedTest
    @MethodSource("test.integration.cli.CliScenario#scenarios")
    public void zip_file_contains_all_files(String name) throws Exception {
        var outputDir = new File(CliScenario.TEST_OUT_DIR).toPath();
        if(!Files.exists(outputDir)) {
            Files.createDirectories(outputDir);
        }

        var scenario = new CliScenario(name, new String[] {
                "fastqc.output_dir=" + outputDir
        });
        
        // execute the scenario, checking its output and that it is successful
        Cli.Execute(scenario).assertSuccess();

        // locate the zip file and extract it
        assertFileExists(scenario.ZipFilePath);
        TestHelpers.unzip(scenario.ZipFilePath, scenario.ZipExtractionPath);

        // locate the inner folder and check it exists
        var folder = new File(scenario.ZipExtractionPath, scenario.ZipInnerFolderName).getAbsolutePath().toString();
        assertFolderExists(folder);

        // assert files
        assertFileExists(folder + "/fastqc_data.txt");
        assertFileExists(folder + "/fastqc_report.html");
        assertFileExists(folder + "/fastqc.fo");
        assertFileExists(folder + "/summary.txt");

        // assert images
        final String ImagesFolderName = "Images";
        var imagesFolder = folder + "/" + ImagesFolderName;
        assertFolderExists(imagesFolder);
        assertImageExists(imagesFolder + "/adapter_content");
        assertImageExists(imagesFolder + "/duplication_levels");
        assertImageExists(imagesFolder + "/per_base_n_content");
        assertImageExists(imagesFolder + "/per_base_quality");
        assertImageExists(imagesFolder + "/per_base_sequence_content");
        assertImageExists(imagesFolder + "/per_sequence_gc_content");
        assertImageExists(imagesFolder + "/per_sequence_quality");
        assertImageExists(imagesFolder + "/sequence_length_distribution");

        // assert icons
        final String IconsFolderName = "Icons";
        var iconsFolder = folder + "/" + IconsFolderName;
        assertFolderExists(iconsFolder);
        assertFileExists(iconsFolder + "/error.png");
        assertFileExists(iconsFolder + "/fastqc_icon.png");
        assertFileExists(iconsFolder + "/tick.png");
        assertFileExists(iconsFolder + "/warning.png");
    }

    private void assertImageExists(String path) {
        assertFileExists(path + ".png");
        assertFileExists(path + ".svg");
    }

    private void assertIconExists(String path) {
        assertFileExists(path + ".png");
    }
}
