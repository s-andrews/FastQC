package test.ui.Modules;

import java.io.File;
import java.io.IOException;
import java.awt.*;
import java.awt.image.BufferedImage;

import org.junit.jupiter.api.Test;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import uk.ac.babraham.FastQC.Modules.BasicStats;
import uk.ac.babraham.FastQC.Sequence.SequenceFactory;
import org.approvaltests.Approvals;
import org.approvaltests.awt.AwtApprovals;
import org.approvaltests.core.Options;
import test.data.TestCases;

public class BasicStatsTest {

    @ParameterizedTest
    @MethodSource("test.data.TestCases#all")
    public void testRendering(String name) throws Exception {
        // arrange
        var fastqFile = TestCases.getTestFastQFile(name);
        var sequenceFile = SequenceFactory.getSequenceFile(fastqFile);
        var stats = new BasicStats();
        while (sequenceFile.hasNext()) {
            var sequence = sequenceFile.next();
            stats.processSequence(sequence);
        }
        var panel = stats.getResultsPanel();
        panel.setSize(200, 200);

        // act
        var image = new BufferedImage(panel.getWidth(), panel.getHeight(), BufferedImage.TYPE_INT_ARGB);
        var g2d = image.createGraphics();
        panel.printAll(g2d);
        g2d.dispose();

        // assert
        AwtApprovals.verify(panel, Approvals.NAMES.withParameters(name));
    }
}