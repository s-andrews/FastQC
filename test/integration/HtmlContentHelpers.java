package test.integration;

import java.util.regex.Pattern;

import org.jsoup.Jsoup;
import org.jsoup.nodes.Document;

/**
 * Utilities for scrubbing dates and normalizing HTML content.
 */
public final class HtmlContentHelpers {
    private HtmlContentHelpers() {}

    // Matches: Mon 22 Sep 2025 (or Tue/Wed..., Jan-Dec)
    private static final Pattern DATE_PATTERN = Pattern.compile(
            "(?:Mon|Tue|Wed|Thu|Fri|Sat|Sun)\\s+\\d{1,2}\\s+(?:Jan|Feb|Mar|Apr|May|Jun|Jul|Aug|Sep|Oct|Nov|Dec)\\s+\\d{4}");

    private static final String PLACEHOLDER = "Mon 01 Jan 2000";

    /** Replace all occurrences of the date pattern with <DATE>. */
    public static String scrubDates(String input) {
        if (input == null)
            return null;
        return DATE_PATTERN.matcher(input).replaceAll(PLACEHOLDER);
    }

    public static String normalizeHtml(String html) {
        Document doc = Jsoup.parse(html);
        doc.outputSettings()
                .prettyPrint(true)
                .outline(false)
                .syntax(Document.OutputSettings.Syntax.html);
        return doc.outerHtml();
    }
}