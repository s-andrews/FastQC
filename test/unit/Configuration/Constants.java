package test.unit.Configuration;

import java.util.Set;

public final class Constants {
    private Constants() {
    }
    public static final String DUPLICATION = "duplication";
    public static final String KMER = "kmer";
    public static final String N_CONTENT = "n_content";
    public static final String OVERREPRESENTED = "overrepresented";
    public static final String QUALITY_BASE = "quality_base";
    public static final String SEQUENCE = "sequence";
    public static final String GC_SEQUENCE = "gc_sequence";
    public static final String QUALITY_SEQUENCE = "quality_sequence";
    public static final String TILE = "tile";
    public static final String SEQUENCE_LENGTH = "sequence_length";
    public static final String ADAPTER = "adapter";

    public static final Set<String> MODULE_NAMES = Set.of(
            KMER,
            N_CONTENT,
            OVERREPRESENTED,
            QUALITY_BASE,
            SEQUENCE,
            GC_SEQUENCE,
            QUALITY_SEQUENCE,
            TILE,
            SEQUENCE_LENGTH,
            ADAPTER);
}