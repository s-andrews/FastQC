package test.integration.models;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Enumeration;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class TestZipFileModel {
    private final File file;
    private final TestScenario scenario;

    public TestZipFileModel(TestScenario scenario, File file) {
        this.scenario = scenario;
        this.file = file;        
    }

    public TestZipFolderModel unzip() throws IOException {
        var targetDir = new File(scenario.ZipExtractionPath);
        if (!targetDir.exists()) {
            Files.createDirectories(targetDir.toPath());
        }

        try (ZipFile zipFile = new ZipFile(file)) {
            Enumeration<? extends ZipEntry> entries = zipFile.entries();

            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                File outFile = new File(targetDir, entry.getName());

                if (entry.isDirectory()) {
                    outFile.mkdirs();
                } else {
                    // Ensure parent directories exist
                    outFile.getParentFile().mkdirs();

                    try (var in = zipFile.getInputStream(entry); var out = new FileOutputStream(outFile)) {
                        in.transferTo(out);
                    }
                }
            }
        }
        var innerFolder = new File(targetDir, scenario.Name + "_fastqc");
        
        return new TestZipFolderModel(scenario, innerFolder);
    }
}
