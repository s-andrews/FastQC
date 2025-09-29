package test.unit.Modules;

import org.junit.jupiter.api.Test;
import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import uk.ac.babraham.FastQC.Modules.BasicStats;

public class BasicStatsTest {

    @Test
    public void testFormatLength_dot() {
        assertEquals("0 bp",    BasicStats.formatLength(0));
        // Should not keep trailing .0
        assertEquals("1 kbp",   BasicStats.formatLength(1_000));
        assertEquals("1 Mbp",   BasicStats.formatLength(1_000_000));
        assertEquals("1 Gbp",   BasicStats.formatLength(1_000_000_000));
        // Should keep one decimal place
        assertEquals("1 bp",        BasicStats.formatLength(1));
        assertEquals("11 bp",       BasicStats.formatLength(11));
        assertEquals("111 bp",      BasicStats.formatLength(111));
        assertEquals("1.1 kbp",     BasicStats.formatLength(1_111));
        assertEquals("11.1 kbp",    BasicStats.formatLength(11_111));
        assertEquals("111.1 kbp",   BasicStats.formatLength(111_111));
        assertEquals("1.1 Mbp",     BasicStats.formatLength(1_111_111));
        assertEquals("11.1 Mbp",    BasicStats.formatLength(11_111_111));
        assertEquals("111.1 Mbp",   BasicStats.formatLength(111_111_111));
        assertEquals("1.1 Gbp",     BasicStats.formatLength(1_111_111_111));
        assertEquals("11.1 Gbp",    BasicStats.formatLength(11_111_111_111L));
        assertEquals("111.1 Gbp",   BasicStats.formatLength(111_111_111_111L));
    }
}