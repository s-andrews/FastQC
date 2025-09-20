package test.integration.models;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertTrue;

public class ExecutionResult {
    private final TestScenario scenario;
    private final int exitCode;
    private final String output;

    public ExecutionResult(TestScenario scenario, int exitCode, String output) {
        this.scenario = scenario;
        this.exitCode = exitCode;
        this.output = output;
    }

    public ExecutionResult AssertExitCode(int expected) {
        assertEquals(expected, exitCode, () -> "Unexpected exit code. Expected " + expected + " but got " + exitCode);
        return this;
    }

    public ExecutionResult AssertExitCodeIsZero() {
        return AssertExitCode(0);
    }

    public ExecutionResult AssertStarted() {
        return AssertOutputContains("Started analysis of " + scenario.FastqFileName);
    }

    public ExecutionResult AssertCompleted() {
        return AssertOutputContains("Analysis complete for " + scenario.FastqFileName);
    }

    public ExecutionResult AssertOutputContains(String expected) {
        assertTrue(output.contains(expected), () -> "Expected output to contain:\n" + expected + "\nGot:\n" + output);
        return this;
    }
}
