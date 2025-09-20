package test.integration;

import org.junit.jupiter.api.Test;

import test.integration.models.ExecutionHelper;
import test.integration.models.TestScenario;

public class OutDirTest {

    @Test
    public void shows_error_when_0() throws Exception {
        ExecutionHelper
                .Execute(new TestScenario(null, new String[] { "fastqc.output_dir=this_folder_does_not_exist" }))
                .AssertExitCode(1)
                .AssertOutputContains("Output dir this_folder_does_not_exist doesn't exist or isn't writeable");
    }
}
