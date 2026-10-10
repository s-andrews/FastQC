package test.unit.Modules;

import java.io.BufferedReader;
import java.io.File;
import java.io.InputStreamReader;
import java.util.Arrays;
import java.util.concurrent.TimeUnit;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

import uk.ac.babraham.FastQC.Modules.GCModel.GCModel;
import uk.ac.babraham.FastQC.Modules.GCModel.GCModelValue;
import uk.ac.babraham.FastQC.Modules.PerSequenceGCContent;
import uk.ac.babraham.FastQC.Sequence.Sequence;

public class GCModelTest {
    @Test
    public void preservesBoundaryAndInteriorWeights() {
        GCModel model = new GCModel(1000);
        assertValue(model.getModelValues(0), 0, 1d / 6);
        assertValue(model.getModelValues(500), 50, 1d / 11);
        assertValue(model.getModelValues(1000), 100, 1d / 6);

        GCModelValue[] boundary = new GCModel(100).getModelValues(50);
        assertEquals(2, boundary.length);
        assertEquals(50, boundary[0].percentage());
        assertEquals(51, boundary[1].percentage());
        assertEquals(0.5, boundary[0].increment());
        assertEquals(0.5, boundary[1].increment());
    }

    @Test
    public void normalizesEveryPercentageAcrossPossibleCounts() {
        for (int length : new int[] {1, 2, 10, 99, 100, 101, 1000, 10000}) {
            GCModel model = new GCModel(length);
            double[] totals = new double[101];
            for (int count = 0; count <= length; count++) {
                for (GCModelValue value : model.getModelValues(count)) {
                    totals[value.percentage()] += value.increment();
                }
            }
            for (double total : totals) assertEquals(1d, total, 1e-12);
        }
    }

    @Test
    public void keepsCountArraysIndependentlyMutable() {
        GCModel model = new GCModel(1000);
        GCModelValue[] first = model.getModelValues(500);
        GCModelValue[] second = model.getModelValues(501);
        assertNotSame(first, second);
        first[0] = new GCModelValue(99, 0.25);
        assertSame(first, model.getModelValues(500));
        assertValue(second, 50, 1d / 11);
    }

    @Test
    public void handlesDiverseLongReadLengthsWithin192MiB() throws Exception {
        // Run separately so the heap bound is independent of the test runner.
        String childClasspath = new File(GCModelTest.class.getProtectionDomain().getCodeSource().getLocation().toURI()).getPath()
            + File.pathSeparator + new File(GCModel.class.getProtectionDomain().getCodeSource().getLocation().toURI()).getPath()
            + File.pathSeparator + System.getProperty("java.class.path");
        Process child = new ProcessBuilder(
            new File(System.getProperty("java.home"), "bin/java").getPath(),
            "-Xmx192m", "-cp", childClasspath,
            LongReadProbe.class.getName()).redirectErrorStream(true).start();
        boolean finished = child.waitFor(30, TimeUnit.SECONDS);
        if (!finished) child.destroyForcibly();
        assertTrue(finished, "Long-read GC analysis timed out");
        StringBuilder output = new StringBuilder();
        try (BufferedReader reader = new BufferedReader(new InputStreamReader(child.getInputStream()))) {
            String line;
            while ((line = reader.readLine()) != null) output.append(line).append('\n');
        }
        assertEquals(0, child.exitValue(), output.toString());
        assertTrue(output.toString().contains("Processed 100 reads"), output.toString());
    }

    private static void assertValue(GCModelValue[] values, int percentage, double increment) {
        assertEquals(1, values.length);
        assertEquals(percentage, values[0].percentage());
        assertEquals(increment, values[0].increment());
    }

    public static class LongReadProbe {
        public static void main(String[] args) {
            PerSequenceGCContent module = new PerSequenceGCContent();
            for (int length = 1000; length <= 100000; length += 1000) {
                char[] bases = new char[length];
                Arrays.fill(bases, 'A');
                module.processSequence(new Sequence(null, new String(bases), null, "read" + length));
            }
            // Keep the module reachable through the end of the probe.
            System.out.println("Processed 100 reads: " + module.name());
        }
    }
}
