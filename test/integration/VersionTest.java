package test.integration;

import org.junit.jupiter.api.Test;

import test.integration.cli.Cli;
import test.integration.cli.CliScenario;

public class VersionTest {

    @Test
    public void displays_name_and_semver() throws Exception {
        Cli
            .Execute(new CliScenario(null, new String[] {"fastqc.show_version=true"}))
            .assertSuccess()
            .assertOutputContains("FastQC v0.12.1");
    }

}
