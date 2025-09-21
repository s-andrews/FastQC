package test.integration;

import java.io.File;
import java.nio.file.Files;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import test.integration.cli.Cli;
import test.integration.cli.CliScenario;

public class OutputTest {

    @ParameterizedTest
    @MethodSource("test.integration.cli.CliScenario#scenarios")
    public void shows_progress_when_processing(String name) throws Exception {
        var outputDir = new File(CliScenario.TEST_OUT_DIR).toPath();
        if (!Files.exists(outputDir)) {
            Files.createDirectories(outputDir);
        }

        var scenario = new CliScenario(name, new String[] {
            "fastqc.output_dir=" + outputDir
        });

        Cli.Execute(scenario)
            .assertOutputContains("Started analysis of " + scenario.FastqFileName)
            .assertOutputContains("Analysis complete for " + scenario.FastqFileName)
            .assertSuccess();
    }
}
