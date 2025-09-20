package test.integration.models;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;

import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.*;

public class TestZipFolderModel {
    private final File folder;
    private final TestScenario scenario;

    private final String ImagesFolderName = "Images";
    private final String IconsFolderName = "Icons";

    public TestZipFolderModel(TestScenario scenario, File folder) {
        this.scenario = scenario;
        this.folder = folder;
    }

    public TestZipFolderModel assertZipFileContent() {
        return assertFileExists("fastqc_data.txt")
                .assertFileExists("fastqc_report.html")
                .assertFileExists("fastqc.fo")
                .assertFileExists("summary.txt")              
                .assertFolderExists(ImagesFolderName)
                .assertFolderExists(IconsFolderName)              
                .assertImageExists("adapter_content")
                .assertImageExists("duplication_levels")
                .assertImageExists("per_base_n_content")
                .assertImageExists("per_base_quality")
                .assertImageExists("per_base_sequence_content")
                .assertImageExists("per_sequence_gc_content")
                .assertImageExists("per_sequence_quality")
                .assertImageExists("sequence_length_distribution");
    }

    public TestZipFolderModel assertFastQcFileMatches() {
        var file = new File(folder, "fastqc_data.txt");
        try {
            var content = Files.readString(file.toPath());
            var comparison = scenario.GetFastQcDataSnapshotFileContent();
            assertEquals(comparison, content, () -> "fastqc_data.txt content does not match snapshot");
        } catch (IOException e) {
            throw new RuntimeException("Failed to read fastqc_data.txt content", e);
        }
        return this;
    }

    public TestZipFolderModel assertFolderExists(String name) {
        var folder = new File(this.folder, name);
        assertTrue(folder.exists(), () -> "Expected folder to exist: " + folder.getAbsolutePath());
        assertTrue(folder.isDirectory(), () -> "Expected path to be a folder: " + folder.getAbsolutePath());
        return this;
    }

    private TestZipFolderModel assertImageExists(String name) {
        return assertFileExists(ImagesFolderName + "/" + name + ".png")
                .assertFileExists(ImagesFolderName + "/" + name + ".svg");
    }

    private TestZipFolderModel assertFileExists(String name) {
        var file = new File(folder, name);
        assertTrue(file.exists(), () -> "Expected file to exist: " + file.getAbsolutePath());
        assertTrue(file.isFile(), () -> "Expected path to be a file: " + file.getAbsolutePath());
        return this;
    }
}