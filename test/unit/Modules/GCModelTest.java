package test.unit.Modules;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.*;

import uk.ac.babraham.FastQC.Modules.GCModel.GCModel;
import uk.ac.babraham.FastQC.Modules.GCModel.GCModelValue;

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

    private static void assertValue(GCModelValue[] values, int percentage, double increment) {
        assertEquals(1, values.length);
        assertEquals(percentage, values[0].percentage());
        assertEquals(increment, values[0].increment());
    }
}
