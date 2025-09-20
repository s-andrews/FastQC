package test.integration;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;
import static org.junit.jupiter.api.Assertions.fail;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import test.integration.models.ExecutionHelper;
import test.integration.models.TestScenario;

public class ThreadsTest {

    @Test
    public void shows_error_when_0() throws Exception {
        ExecutionHelper
            .Execute(new TestScenario(null, new String[] {"fastqc.threads=0"}))
            .AssertExitCode(1)
            .AssertOutputContains("Number of threads must be >= 1");
    }

    @Test
    public void shows_error_when_negative() throws Exception {
        ExecutionHelper
            .Execute(new TestScenario(null, new String[] {"fastqc.threads=-1"}))
            .AssertExitCode(1)
            .AssertOutputContains("Number of threads must be >= 1");
    }
    
    // TODO: 
    // The test is disabled for now as invalid parameter value currently results in an exception and stack-trace
    // Ideally, it should show a user-friendly error message instead
    // It fails because the input value is assumed to be an integer
    // threads = Integer.parseInt(System.getProperty("fastqc.threads"));
    @Test
    @Disabled
    public void shows_error_when_alpha() throws Exception {
        ExecutionHelper
            .Execute(new TestScenario(null, new String[] {"fastqc.threads=abc"}))
            .AssertExitCode(1)
            .AssertOutputContains("Number of threads must be >= 1");
    }
}
