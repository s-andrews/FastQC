package test.unit.Modules;

import static org.junit.jupiter.api.Assertions.assertArrayEquals;

import java.lang.reflect.Field;
import java.lang.reflect.Method;
import java.util.Random;

import org.junit.jupiter.api.Test;

import uk.ac.babraham.FastQC.Modules.AdapterContent;

/**
 * The single pass adapter search must give exactly the same first match
 * positions as String.indexOf for every adapter.
 */
public class AdapterContentTest {

    private static final String [] ADAPTERS = {
        "AGATCGGAAGAG", "TGGAATTCTCGG", "GATCGTCGGACT", "CTGTCTCTTATA", "AAAAAAAAAAAA", "GGGGGGGGGGGG"
    };

    private final AdapterContent module = new AdapterContent();

    private int [] findAdapters (String seq) throws Exception {
        Method find = AdapterContent.class.getDeclaredMethod("findAdapters", String.class);
        find.setAccessible(true);
        find.invoke(module, seq);
        Field hits = AdapterContent.class.getDeclaredField("firstHit");
        hits.setAccessible(true);
        return ((int []) hits.get(module)).clone();
    }

    private static int [] indexOf (String seq) {
        int [] expected = new int[ADAPTERS.length];
        for (int a=0;a<ADAPTERS.length;a++) expected[a] = seq.indexOf(ADAPTERS[a]);
        return expected;
    }

    private void check (String seq) throws Exception {
        assertArrayEquals(indexOf(seq), findAdapters(seq), seq);
    }

    @Test
    public void testEdgeCases () throws Exception {
        check("");
        check("AGATCGGAAGA");
        check("AGATCGGAAGAG");
        check("NAGATCGGAAGAGN");
        check("AGATCGGANAGAG");
        check("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA");
        check("AAAAAAAAAAANAAAAAAAAAAAAGGGGGGGGGGG.GGGGGGGGGGGGG");
        check("CTGTCTCTTATAAGATCGGAAGAGTGGAATTCTCGGGATCGTCGGACT");
        check("agatcggaagag");
        check("AGATCGGAAGAŇ");
    }

    @Test
    public void testRandomReads () throws Exception {
        Random random = new Random(42);
        String bases = "ACGTACGTACGTACGTN";
        for (int r=0;r<100_000;r++) {
            StringBuilder sb = new StringBuilder();
            int length = random.nextInt(160);
            while (sb.length() < length) {
                if (random.nextInt(20) == 0) {
                    // plant a whole or partial adapter
                    String adapter = ADAPTERS[random.nextInt(ADAPTERS.length)];
                    sb.append(adapter, 0, 1+random.nextInt(adapter.length()));
                }
                else {
                    sb.append(bases.charAt(random.nextInt(bases.length())));
                }
            }
            check(sb.toString());
        }
    }
}
