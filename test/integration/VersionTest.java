package test.integration;

import org.junit.jupiter.api.Test;

import test.integration.models.ExecutionHelper;
import test.integration.models.TestScenario;

public class VersionTest {

    @Test
    public void show_version_displays_name_and_semver() throws Exception {
        ExecutionHelper
                .Execute(new TestScenario(null, new String[] {"fastqc.show_version=true"}))
                .AssertExitCodeIsZero()
                .AssertOutputContains("FastQC v0.12.1");
    }

}
