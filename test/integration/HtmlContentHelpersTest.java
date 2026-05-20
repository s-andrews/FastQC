package test.integration;

import static org.junit.jupiter.api.Assertions.assertEquals;

import org.junit.jupiter.api.Test;

public class HtmlContentHelpersTest {

    @Test
    public void matches_and_replaces_date_substring() {
        String input = "<div id=\"header_filename\">Fri 25 Dec 2025<br>complex.fastq</div>";
        String expected = "<div id=\"header_filename\">Mon 01 Jan 2000<br>complex.fastq</div>";
        String actual = HtmlContentHelpers.scrubDates(input);
        assertEquals(expected, actual);
    }
}
