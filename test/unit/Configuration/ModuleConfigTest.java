package test.unit.Configuration;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static test.unit.Configuration.Constants.*;

import java.util.Set;

import uk.ac.babraham.FastQC.Modules.ModuleConfig;

public class ModuleConfigTest {

    private ModuleConfig config;

    @BeforeEach
    public void setUp() {
        config = new ModuleConfig();
    }

    @Test
    public void testEnabledModules() {
        for (var moduleName : Set.of(
                N_CONTENT,
                OVERREPRESENTED,
                QUALITY_BASE,
                SEQUENCE,
                GC_SEQUENCE,
                QUALITY_SEQUENCE,
                TILE,
                SEQUENCE_LENGTH,
                ADAPTER)) {
            assertModuleEnabled(moduleName);
        }
    }

    @Test
    public void testDisabledModules() {
        assertModuleDisabled(KMER);
    }

    private void assertModuleEnabled(String moduleName) {
        assertEquals(0, config.getParam(moduleName, "ignore"),
                "Module " + moduleName + " should be enabled by default");
    }

    private void assertModuleDisabled(String moduleName) {
        assertEquals(1, config.getParam(moduleName, "ignore"),
                "Module " + moduleName + " should be disabled by default");
    }

}