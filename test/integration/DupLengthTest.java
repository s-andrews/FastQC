package test.integration;

import java.io.File;
import java.nio.file.Files;
import java.nio.file.Path;

import org.junit.jupiter.api.Test;

import test.integration.cli.Cli;
import test.integration.cli.CliScenario;

public class DupLengthTest {

    @Test
    public void succeeds_when_dup_length_is_longer_than_reads() throws Exception {
        Path outputDir = new File(CliScenario.TEST_OUT_DIR).toPath();
        if (!Files.exists(outputDir)) {
            Files.createDirectories(outputDir);
        }

        CliScenario scenario = new CliScenario("minimal", new String[] {
                "fastqc.output_dir=" + outputDir,
                "fastqc.dup_length=20"
        });

        Cli.Execute(scenario).assertSuccess();
    }
}
