package test.unit;

import java.awt.Font;
import java.awt.GraphicsEnvironment;
import javax.swing.UIManager;
import org.junit.jupiter.api.Test;

public class FontTests {

    @Test
    public void listAllAvailableFonts() {
        System.out.println("=== ALL AVAILABLE FONTS ===");
        for (Font font : GraphicsEnvironment.getLocalGraphicsEnvironment().getAllFonts()) {
            System.out.println("Font: " + font.getName() + " (Family: " + font.getFamily() + ")");
        }
    }

    @Test
    public void showDefaultFonts() {
        System.out.println("\n=== DEFAULT SWING FONTS ===");
        
        // Show key Swing font settings
        System.out.println("Default Button Font: " + UIManager.getFont("Button.font"));
        System.out.println("Default Label Font: " + UIManager.getFont("Label.font"));
        System.out.println("Default Panel Font: " + UIManager.getFont("Panel.font"));
        System.out.println("Default TextField Font: " + UIManager.getFont("TextField.font"));
        
        // Test logical font mappings
        System.out.println("\n=== LOGICAL FONT MAPPINGS ===");
        Font sansSerif = new Font(Font.SANS_SERIF, Font.PLAIN, 12);
        Font serif = new Font(Font.SERIF, Font.PLAIN, 12);
        Font monospaced = new Font(Font.MONOSPACED, Font.PLAIN, 12);
        Font dialog = new Font(Font.DIALOG, Font.PLAIN, 12);
        
        System.out.println("Font.SANS_SERIF maps to: " + sansSerif.getName() + " (Family: " + sansSerif.getFamily() + ")");
        System.out.println("Font.SERIF maps to: " + serif.getName() + " (Family: " + serif.getFamily() + ")");
        System.out.println("Font.MONOSPACED maps to: " + monospaced.getName() + " (Family: " + monospaced.getFamily() + ")");
        System.out.println("Font.DIALOG maps to: " + dialog.getName() + " (Family: " + dialog.getFamily() + ")");
    }
}
