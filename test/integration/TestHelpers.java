package test.integration;

import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.nio.file.Files;
import java.util.Enumeration;
import java.util.regex.Pattern;
import java.util.zip.ZipEntry;
import java.util.zip.ZipFile;

public class TestHelpers {
    
    public static void assertFolderExists(String path) {
        File folder = new File(path);
        assertTrue(folder.exists(), () -> "Expected folder to exist: " + folder.getAbsolutePath());
        assertTrue(folder.isDirectory(), () -> "Expected path to be a folder: " + folder.getAbsolutePath());
    }

    public static void assertFileExists(String path) {
        File file = new File(path);
        assertTrue(file.exists(), () -> "Expected file to exist: " + file.getAbsolutePath());
        assertTrue(file.isFile(), () -> "Expected path to be a file: " + file.getAbsolutePath());
    }
    
    public static void unzip(String zipFilePath, String extractionPath) throws IOException {
        File targetDir = new File(extractionPath);
        if (!targetDir.exists()) {
            Files.createDirectories(targetDir.toPath());
        }

        try (ZipFile zipFile = new ZipFile(zipFilePath)) {
            Enumeration<? extends ZipEntry> entries = zipFile.entries();

            while (entries.hasMoreElements()) {
                ZipEntry entry = entries.nextElement();
                File outFile = new File(targetDir, entry.getName());

                if (entry.isDirectory()) {
                    outFile.mkdirs();
                } else {
                    // Ensure parent directories exist
                    outFile.getParentFile().mkdirs();

                    try (InputStream in = zipFile.getInputStream(entry); FileOutputStream out = new FileOutputStream(outFile)) {
                        in.transferTo(out);
                    }
                }
            }
        }
    }
    
    private static final Pattern DATE_PATTERN = Pattern.compile("\\b\\w{3}\\s+\\d{1,2}\\s+\\w{3}\\s+\\d{4}\\b");

    public static String removeDates(String input) {
        return DATE_PATTERN.matcher(input).replaceAll("").trim();
    }
}
