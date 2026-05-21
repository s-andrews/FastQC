package test.integration;

import org.junit.jupiter.api.Test;

import test.integration.cli.Cli;
import test.integration.cli.CliScenario;
import uk.ac.babraham.FastQC.FastQCApplication;

public class VersionTest {

    @Test
    public void displays_name_and_semver() throws Exception {
        Cli
            .Execute(new CliScenario(null, new String[] {"fastqc.show_version=true"}))
            .assertSuccess()
            .assertOutputContains("FastQC v" + FastQCApplication.VERSION);
    }

}
