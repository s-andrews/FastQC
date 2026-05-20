# Tests
- [Test Data](#test-data)
  - [Use with Parameterised Tests](#use-with-parameterised-tests)
- [Running Tests](#running-tests)
  - [Unit Tests](#unit-tests)
  - [Integration Tests](#integration-tests)
  - [UI Tests](#ui-tests)

## Test Data

Test data is located here `/test/data`

Test scenarios are limited at the moment, just two hand crafted fastq files
- `simple` a single read, 16bp long.
- `complex` a collection of five reads, 16bp each, 80bp in total.

### Use with Parameterised Tests
The test data is designed to be used with parametised tests.
It consists of two collections at the moment, `fast` and `slow` as well as the aggregate of those called `all`. Fast tests are handy to use interactively while writing tests or diagnosing issues. The slow tests will operate on real world data and take much longer to run and will more likely form part of a CI / regression testing process.

```java
    public static Stream<String> all() {
        // return the union of fast and slow test cases
        return Stream.concat(fast(), slow());
    }
    public static Stream<String> fast() {
        return Stream.of("minimal");
    }
    public static Stream<String> slow() {
        return Stream.of("complex");
    }
```
These can be consumed in your test per the below:

```java
public class YourTest {

    @ParameterizedTest
    @MethodSource("test.data.TestCases#all")
    public void testSomething(String name) throws Exception {
        var fastqFile = TestCases.getTestFastQFile(name);
        ...
```

## Running Tests

### Unit Tests
```sh
ant clean build unit-test
```

### Integration Tests
```sh
ant clean build integration-test
```

### UI Tests

Running ui tests is slightly more complicated, but the dev container installs the dependencies required to simplify this.
```sh
xvfb-run -a ant clean build ui-test
```