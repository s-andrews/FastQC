package test.integration.models;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Parameter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.stream.Stream;

public class TestScenario {

    private final String TEST_DATA_DIR = "test/integration/data/";
    public static final String TEST_OUT_DIR = "test/integration/output/";

    public TestScenario(String[] parameters) {
        Parameters = parameters;
    }

    public TestScenario(String name, String[] parameters) {
        Name = name;
        Parameters = parameters;
        FastqFileName = name + ".fastq";
        FastqFilePath = TEST_DATA_DIR + FastqFileName;
        ZipFilePath = TEST_OUT_DIR + name + "_fastqc.zip";
        ZipExtractionPath = TEST_OUT_DIR;
        if(!Files.exists(Path.of(TEST_OUT_DIR))) {
            try {
                Files.createDirectories(Path.of(TEST_OUT_DIR));
            } catch (IOException e) {
                throw new RuntimeException("Failed to create output directory " + TEST_OUT_DIR, e);
            }
        }
    }

    public String Name;
    public String FastqFileName;
    public String FastqFilePath;
    public String ZipFilePath;
    public String ZipExtractionPath;
    public String[] Parameters;


    public TestZipFileModel GetZipFile() {
        var file = new java.io.File(ZipFilePath);
        assertValidFile(file);
        return new TestZipFileModel(this, file);
    }

    public String GetFastQcDataSnapshotFileContent() throws java.io.IOException {
        var path = new java.io.File(TEST_DATA_DIR + Name + "_fastqc_data.txt").toPath();
        try {
            return Files.readString(path);
        } catch (IOException e) {
            throw new RuntimeException("Failed to read file content of " + path.getFileName(), e);
        }
    }

    public void assertValidFile(File file) {
        assertTrue(file.exists(), "Expected file " + file.getPath() + " to exist");
        assertTrue(file.length() > 0, "Expected file " + file.getPath() + " to have size > 0");
    }
}