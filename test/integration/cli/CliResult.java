package test.integration.cli;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertNotEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class CliResult {
    private final CliScenario scenario;
    private final int exitCode;
    private final String output;
    private final int successfulExitCode = 0;

    public CliResult(CliScenario scenario, int exitCode, String output) {
        this.scenario = scenario;
        this.exitCode = exitCode;
        this.output = output;
    }

    public CliResult assertSuccess() {
        assertEquals(successfulExitCode, exitCode, () -> "Unexpected exit code. Expected " + successfulExitCode + " but got " + exitCode);
        return this;
    }

    public CliResult assertFailure() {
        assertNotEquals(successfulExitCode, exitCode, () -> "Unexpected exit code. Expected NOT " + successfulExitCode + " but got " + exitCode);
        return this;
    }

    public CliResult assertOutputContains(String expected) {
        assertTrue(output.contains(expected), () -> "Expected output to contain:\n" + expected + "\nGot:\n" + output);
        return this;
    }
}
