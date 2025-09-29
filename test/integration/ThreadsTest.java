package test.integration;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.api.Test;

import test.integration.cli.Cli;
import test.integration.cli.CliScenario;

public class ThreadsTest {

    @Test
    public void shows_error_when_value_is_zero() throws Exception {
        Cli
            .Execute(new CliScenario(new String[] {"fastqc.threads=0"}))
            .assertFailure()
            .assertOutputContains("Number of threads must be >= 1");
    }

    @Test
    public void shows_error_when_value_is_negative() throws Exception {
        Cli
            .Execute(new CliScenario(new String[] {"fastqc.threads=-1"}))
            .assertFailure()
            .assertOutputContains("Number of threads must be >= 1");
    }
    
    // TODO: 
    // The test is disabled for now as invalid parameter value currently results in an exception and stack-trace
    // Ideally, it should show a user-friendly error message instead
    // It fails because the input value is assumed to be an integer
    // threads = Integer.parseInt(System.getProperty("fastqc.threads"));
    @Test
    @Disabled
    public void shows_error_when_value_is_non_numeric() throws Exception {
        Cli
            .Execute(new CliScenario(new String[] {"fastqc.threads=abc"}))
            .assertFailure()
            .assertOutputContains("Number of threads must be >= 1");
    }
}
