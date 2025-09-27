package test.data;

import java.io.File;
import java.io.IOException;
import java.util.stream.Stream;

public class TestCases {

    public static String TEST_DATA_DIR = "test/data/";

    /// Returns a File object for the specified test FastQ file
    public static File getTestFastQFile(String testCaseName) throws IOException {
        var file = new File(TEST_DATA_DIR + testCaseName + ".fastq");
        if (!file.exists()) {
            throw new IOException("Test data file for test case " + testCaseName + " not found: " + file.getAbsolutePath());
        }
        return file;
    }

    /// Returns all test case names (without extensions)
    public static Stream<String> all() {
        // return the union of fast and slow test cases
        return Stream.concat(fast(), slow());
    }

    /// Returns only fast test cases, suitable for interactive testing
    /// (i.e. small files that run quickly)
    public static Stream<String> fast() {
        return Stream.of("minimal");
    }

    /// Returns only slow test cases, suitable for automated testing
    /// (i.e. larger files that take longer to run)
    public static Stream<String> slow() {
        return Stream.of("complex");
    }
}
