package test.integration;

import static org.junit.jupiter.api.Assertions.assertFalse;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;
import org.junit.jupiter.api.Test;

import test.integration.cli.Cli;
import test.integration.cli.CliScenario;

public class TruncatedFileTest {

    // truncated.fastq ends part way through the quality line of its second entry
    @Test
    public void file_cut_off_mid_entry_fails_without_a_report() throws Exception {
        Path outputDir = new File(CliScenario.TEST_OUT_DIR).toPath();
        if (!Files.exists(outputDir)) {
            Files.createDirectories(outputDir);
        }

        CliScenario scenario = new CliScenario("truncated", new String[] {
                "fastqc.output_dir=" + outputDir
        });

        Cli.Execute(scenario)
                .assertFailure()
                .assertOutputContains("Quality string length (8) doesn't match sequence length (16)");

        assertFalse(new File(scenario.ZipFilePath).exists(), "No report should be written for a truncated file");
    }
}
