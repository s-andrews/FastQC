package uk.ac.babraham.FastQC.FileFilters;

import org.junit.Test;
import uk.ac.babraham.FastQC.FileFilters.SequenceFileFilter;
import static org.junit.Assert.*;
import java.io.File;


public class SequenceFileFilterTest {

    private final SequenceFileFilter filter = new SequenceFileFilter();

    @Test
    public void testGetDescription() {
        assertEquals("Sequence Files", filter.getDescription());
    }

    @Test
    public void testAcceptsDirectory() {
        var dir = new File("fake-directory-name") {
            @Override
            public boolean isDirectory() {
                return true;
            }
        };
        assertTrue(filter.accept(dir));
    }

    @Test
    public void testAcceptsValidExtensions() {
        var fileNames = new String[]{
            "file.txt.gz", "file.fastq.gz", "file.fq.gz", "file.fq",
            "file.txt.bz2", "file.fastq.bz2", "file.txt", "file.fastq",
            "file.bam", "file.sam", "file.compact-reads", "file.goby"
        };
        for (String name : fileNames) {
            var f = new File(name);
            assertTrue("Should accept: " + name, filter.accept(f));
        }
    }

    @Test
    public void testRejectsInvalidExtensions() {
        var fileNames = new String[]{
            "file.doc", "file.jpg", "file.fasta", "file.zip", "file.gz", "file.bz2", "file"
        };
        for (String name : fileNames) {
            var f = new File(name);
            assertFalse("Should reject: " + name, filter.accept(f));
        }
    }

    @Test
    public void testCaseInsensitiveExtensions() {
        var f = new File("SAMPLE.FASTQ.GZ");
        assertTrue(filter.accept(f));
    }
}