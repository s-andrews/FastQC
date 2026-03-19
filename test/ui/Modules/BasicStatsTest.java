package test.ui.Modules;

import java.io.File;

import javax.swing.JPanel;

import java.awt.*;
import java.awt.image.BufferedImage;

import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import uk.ac.babraham.FastQC.Modules.BasicStats;
import uk.ac.babraham.FastQC.Sequence.Sequence;
import uk.ac.babraham.FastQC.Sequence.SequenceFactory;
import uk.ac.babraham.FastQC.Sequence.SequenceFile;

import org.approvaltests.Approvals;
import org.approvaltests.awt.AwtApprovals;
import test.data.TestCases;

public class BasicStatsTest {

    @ParameterizedTest
    @MethodSource("test.data.TestCases#all")
    public void testRendering(String name) throws Exception {
        // arrange
        File fastqFile = TestCases.getTestFastQFile(name);
        SequenceFile sequenceFile = SequenceFactory.getSequenceFile(fastqFile);
        BasicStats stats = new BasicStats();
        while (sequenceFile.hasNext()) {
            Sequence sequence = sequenceFile.next();
            stats.processSequence(sequence);
        }
        JPanel panel = stats.getResultsPanel();
        panel.setSize(200, 200);

        // act
        BufferedImage image = new BufferedImage(panel.getWidth(), panel.getHeight(), BufferedImage.TYPE_INT_ARGB);
        Graphics2D g2d = image.createGraphics();
        panel.printAll(g2d);
        g2d.dispose();

        // assert
        AwtApprovals.verify(panel, Approvals.NAMES.withParameters(name));
    }
}