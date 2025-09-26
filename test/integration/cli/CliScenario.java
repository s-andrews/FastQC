package test.integration.cli;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Parameter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.Arrays;
import java.util.stream.Stream;

import test.data.TestCases;
import static test.data.TestCases.TEST_DATA_DIR;

public class CliScenario {

    public static String TEST_OUT_DIR = "test/integration/output/";

    public CliScenario(String[] parameters) {
        Parameters = parameters;
    }

    public CliScenario(String name, String[] parameters) {
        Name = name;
        Parameters = parameters;
        FastqFileName = name + ".fastq";
        FastqFilePath = TEST_DATA_DIR + FastqFileName;
        ZipFilePath = TEST_OUT_DIR + name + "_fastqc.zip";
        ZipExtractionPath = TEST_OUT_DIR;
        ZipInnerFolderName = name + "_fastqc";
        FastqcDataSnapshotPath = TEST_DATA_DIR + name + "_fastqc_data.txt";
        FastqcHtmlSnapshotPath = TEST_DATA_DIR + name + "_fastqc.html";
    }

    public String Name;
    public String FastqFileName;
    public String FastqFilePath;
    public String ZipFilePath;
    public String ZipExtractionPath;
    public String ZipInnerFolderName;
    public String FastqcDataSnapshotPath;
    public String FastqcHtmlSnapshotPath;
    public String[] Parameters;

    @Override
    public String toString() {
        return new StringBuilder()
                .append("Scenario: ").append(Name).append(System.lineSeparator())
                .append("  FastqFile: ").append(FastqFilePath).append(System.lineSeparator())
                .append("  ZipFile:   ").append(ZipFilePath).append(System.lineSeparator())
                .append("  Params:    ").append(String.join(", ", Parameters))
                .toString();
    }
}