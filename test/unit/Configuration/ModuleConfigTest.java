package test.unit.Configuration;

import org.junit.jupiter.api.BeforeEach;
import org.junit.jupiter.api.Test;

import static org.junit.jupiter.api.Assertions.assertEquals;
import static test.unit.Configuration.Constants.*;

import java.util.Set;

import uk.ac.babraham.FastQC.FastQCConfig;
import uk.ac.babraham.FastQC.Modules.ModuleConfig;

public class ModuleConfigTest {

    private ModuleConfig config;

    @BeforeEach
    public void setUp() {
        config = new ModuleConfig(new FastQCConfig());
    }

    @Test
    public void testEnabledModules() {
        // arrange
        var enabledModules = Set.of(
                N_CONTENT,
                OVERREPRESENTED,
                QUALITY_BASE,
                SEQUENCE,
                GC_SEQUENCE,
                QUALITY_SEQUENCE,
                TILE,
                SEQUENCE_LENGTH,
                ADAPTER);
        // act
        for (var moduleName : enabledModules) {
            // assert
            assertModuleEnabled(moduleName);
        }
    }

    @Test
    public void testDisabledModules() {
        // arrange
        var disabledModules = Set.of(KMER);
        // act
        for (var moduleName : disabledModules) {
            // assert
            assertModuleDisabled(moduleName);
        }
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