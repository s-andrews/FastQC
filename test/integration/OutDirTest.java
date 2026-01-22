package test.integration;
import org.junit.jupiter.api.Test;

import test.integration.cli.Cli;
import test.integration.cli.CliScenario;

public class OutDirTest {

    @Test
    public void shows_error_when_output_dir_does_not_exist() throws Exception {
        Cli
            .Execute(new CliScenario(new String[] { "fastqc.output_dir=this_folder_does_not_exist" }))
            .assertFailure()
            .assertOutputContains("Output dir this_folder_does_not_exist doesn't exist or isn't writeable");
    }
}