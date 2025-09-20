package test.integration;

import java.util.stream.Stream;

import org.junit.jupiter.api.Disabled;
import org.junit.jupiter.params.ParameterizedTest;
import org.junit.jupiter.params.provider.MethodSource;

import test.integration.models.ExecutionHelper;
import test.integration.models.TestScenario;

public class ExecutionTest {

        private static final String[] CommonParameters = new String[] {
                        "fastqc.output_dir=" + TestScenario.TEST_OUT_DIR
        };

        public static Stream<TestScenario> scenarios() {
                return Stream.of(
                                new TestScenario("minimal", CommonParameters),
                                new TestScenario("complex", CommonParameters));
        }

        @ParameterizedTest(name = "[{index}] {0}")
        @MethodSource("scenarios")
        public void zip_file_contains_all_files(TestScenario scenario) throws Exception {

                var result = ExecutionHelper.Execute(scenario);

                result.AssertStarted()
                                .AssertCompleted()
                                .AssertExitCodeIsZero();

                var zip = scenario.GetZipFile();
                var folder = zip.unzip();

                folder.assertZipFileContent();
                folder.assertFastQcFileMatches();
        }
}
