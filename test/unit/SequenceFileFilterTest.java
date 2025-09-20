package test.unit;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static org.junit.jupiter.api.Assertions.assertFalse;
import static org.junit.jupiter.api.Assertions.assertTrue;

import java.io.File;

import org.junit.jupiter.api.Test;

import uk.ac.babraham.FastQC.FileFilters.SequenceFileFilter;


public class SequenceFileFilterTest {

    private final SequenceFileFilter filter = new SequenceFileFilter();

    @Test
    public void testGetDescription() {
        assertEquals("Sequence Files", filter.getDescription());
    }

    @Test
    public void testAcceptsDirectory() {
        File dir = new File("fake-directory-name") {
            @Override
            public boolean isDirectory() {
                return true;
            }
        };
        assertTrue(filter.accept(dir));
    }

    @Test
    public void testAcceptsValidExtensions() {
        String[] fileNames = new String[]{
            "file.txt.gz", "file.fastq.gz", "file.fq.gz", "file.fq",
            "file.txt.bz2", "file.fastq.bz2", "file.txt", "file.fastq",
            "file.bam", "file.sam", "file.compact-reads", "file.goby"
        };
        for (String name : fileNames) {
            File f = new File(name);
            assertTrue(filter.accept(f), "Should accept: " + name);
        }
    }

    @Test
    public void testRejectsInvalidExtensions() {
        String[] fileNames = new String[]{
            "file.doc", "file.jpg", "file.fasta", "file.zip", "file.gz", "file.bz2", "file"
        };
        for (String name : fileNames) {
            File f = new File(name);
            assertFalse(filter.accept(f), "Should reject: " + name);
        }
    }

    @Test
    public void testCaseInsensitiveExtensions() {
        File f = new File("SAMPLE.FASTQ.GZ");
        assertTrue(filter.accept(f));
    }
}