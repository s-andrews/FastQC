package test.integration.cli;

import static org.junit.jupiter.api.Assertions.fail;

import java.nio.file.Path;

public class Cli {
    public static final int DEFAULT_TIMEOUT_SECONDS = 30;

    public static CliResult Execute(CliScenario scenario) throws Exception {
        System.out.println("Running scenario: " + scenario);

        // Build a classpath that matches the CLI example plus compiled classes in bin/
        String sep = java.io.File.pathSeparator;
        String cp = String.join(sep,
                "bin",
                ".",
                "jbzip2-0.9.jar",
                "cisd-jhdf5.jar",
                "htsjdk.jar");

        java.util.List<String> cmd = new java.util.ArrayList<>();
        cmd.add("java");
        cmd.add("-Xmx250m");
        for (String parameter : scenario.Parameters) {
            cmd.add("-D" + parameter);
        }
        cmd.add("-cp");
        cmd.add(cp);
        cmd.add("uk.ac.babraham.FastQC.FastQCApplication");
        if (scenario.FastqFilePath != null) {
            cmd.add(scenario.FastqFilePath);
        }
        System.out.println("Command: " + String.join(" ", cmd));

        ProcessBuilder pb = new ProcessBuilder(cmd);
        pb.redirectErrorStream(true);
        // Run from project root so relative classpath entries resolve
        pb.directory(new java.io.File(System.getProperty("user.dir")));

        Process p = pb.start();
        StringBuilder out = new StringBuilder();
        try (java.io.BufferedReader reader = new java.io.BufferedReader(
                new java.io.InputStreamReader(p.getInputStream()))) {
            String line;
            while ((line = reader.readLine()) != null) {
                out.append(line).append(System.lineSeparator());
            }
        }

        boolean finished = p.waitFor(DEFAULT_TIMEOUT_SECONDS, java.util.concurrent.TimeUnit.SECONDS);
        if (!finished) {
            p.destroyForcibly();
            fail("FastQC CLI did not exit within " + DEFAULT_TIMEOUT_SECONDS + "s. Command: " + String.join(" ", cmd));
        }

        var exitCode = p.exitValue();
        var output = out.toString();
        System.out.println("Output:\n" + output); 
        System.out.println("Exit code:\n" + exitCode);
        return new CliResult(scenario, exitCode, output);
    }
}