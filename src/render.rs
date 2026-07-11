//! A module to support rendering of sequences/bases with ANSI escape codes

use std::fmt::Display;

#[derive(Debug, Clone, PartialEq, Eq, Default)]
enum AnsiStyle {
    #[default]
    Normal,
    Faint,
    Italic,
    Bold,
    Underline,
    Strikethrough,
}

impl AnsiStyle {
    pub const fn code(self) -> &'static str {
        match self {
            AnsiStyle::Normal => "0",
            AnsiStyle::Bold => "1",
            AnsiStyle::Faint => "2",
            AnsiStyle::Italic => "3",
            AnsiStyle::Underline => "4",
            AnsiStyle::Strikethrough => "9",
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SeqStyler {
    ansi: bool,
    colour_background: Option<u8>,
    colour_foreground: Option<u8>,
    bold: bool,
    faint: bool,
    underline: bool,
    strikethrough: bool,
    italic: bool,
    fixed_width: bool,
    overwrite: bool,
}

impl Default for SeqStyler {
    fn default() -> Self {
        Self::const_default()
    }
}

impl SeqStyler {
    pub const RESET: &'static str = "\x1b[0m";
    pub const ESC: &'static str = "\x1b";
    pub const CSI: &'static str = "\x1b[";

    const fn const_default() -> Self {
        Self {
            ansi: false,
            colour_background: None,
            colour_foreground: None,
            fixed_width: true,
            overwrite: true,
            bold: false,
            faint: false,
            underline: false,
            strikethrough: false,
            italic: false,
        }
    }

    /// Take any string with a display implementation and
    /// format it according to the SeqStyler settings
    ///
    ///
    /// ```
    /// let s = "ACTGCA"
    ///
    /// // Create Styler
    /// let styler = SeqStyler::new()
    ///     .bold()
    ///     .italic()
    ///     .colour_background(10)
    ///
    /// // Paint style onto any type with a Display impl
    /// let formatted_string = styler::paint(s)
    ///
    /// // Print to any terminal that supports standard ANSI codes and 8bit colour
    /// println!("{formatted_string}")
    ///
    /// ```
    pub fn paint(&self, seq: impl Display) -> String {
        let seq = seq.to_string();

        // Add text padding if fixed_width = TRUE
        let text = if self.fixed_width {
            seq.chars().map(|ch| format!(" {ch} ")).collect::<String>()
        } else {
            seq
        };

        // Return plain text if no ANSI formatting has been requested.
        if !self.ansi {
            return text;
        }

        let codes = self.sgr_codes();

        if codes.is_empty() {
            return text;
        }

        let ansi_prefix = format!("{}{}m", Self::CSI, codes.join(";"));

        format!("{ansi_prefix}{text}{}", Self::RESET)
    }

    // Building up style settings
    pub const fn colour_background(mut self, col: u8) -> Self {
        self.colour_background = Some(col);
        self.ansi = true;
        self
    }

    pub const fn colour_foreground(mut self, col: u8) -> Self {
        self.colour_foreground = Some(col);
        self.ansi = true;
        self
    }

    /// Bold
    pub const fn bold(mut self) -> Self {
        self.bold = true;
        self.ansi = true;
        self
    }

    /// Strikethrough
    pub const fn strikethrough(mut self) -> Self {
        self.strikethrough = true;
        self.ansi = true;
        self
    }

    /// Underline
    pub const fn underline(mut self) -> Self {
        self.underline = true;
        self.ansi = true;
        self
    }

    /// italic
    pub const fn italic(mut self) -> Self {
        self.italic = true;
        self.ansi = true;
        self
    }

    /// faint (dimmed)
    pub const fn faint(mut self) -> Self {
        self.faint = true;
        self.ansi = true;
        self
    }

    /// Instead of overwriting previous styles by sending format reset ANSI string
    /// allow previous styles to stay active
    pub const fn do_not_overwrite(mut self) -> Self {
        self.overwrite = false;
        self.ansi = true;
        self
    }

    /// Disable fixed-width printing
    pub const fn disable_fixed_width(mut self) -> Self {
        self.fixed_width = false;
        self
    }

    // Get all the sgr codes as a vector (we can later semicolon-separate them to format)
    fn sgr_codes(&self) -> Vec<String> {
        let mut codes = Vec::new();

        if self.overwrite {
            codes.push(AnsiStyle::Normal.code().to_string());
        }

        if self.bold {
            codes.push(AnsiStyle::Bold.code().to_string());
        }

        if self.italic {
            codes.push(AnsiStyle::Italic.code().to_string());
        }

        if self.underline {
            codes.push(AnsiStyle::Underline.code().to_string());
        }

        if self.strikethrough {
            codes.push(AnsiStyle::Strikethrough.code().to_string());
        }

        if self.faint {
            codes.push(AnsiStyle::Faint.code().to_string());
        }

        if let Some(col) = self.colour_foreground {
            codes.push(format!("38;5;{col}"));
        }

        if let Some(col) = self.colour_background {
            codes.push(format!("48;5;{col}"));
        }

        codes
    }

    /// Create a new styler
    pub const fn new() -> Self {
        Self::const_default()
    }
}

// Implement Preset Stylers

impl SeqStyler {
    pub const PLAIN: Self = Self::new();

    pub const BORING: Self = Self::new().colour_foreground(15).colour_background(240);

    pub const HIGHLIGHT: Self = Self::new()
        .bold()
        .colour_foreground(16)
        .colour_background(220);

    pub const DIMMED: Self = Self::new()
        .faint()
        .colour_foreground(15)
        .colour_background(240);
}
