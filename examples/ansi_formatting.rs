// This example is less for crate users and more for crate extenders.
//
// Almost everyone reading this wll probably find the 'highlight_sequences' example more useful.
//
//
// There are a lot of different ways we might want to print a sequence to the terminals.
// For example highlighting just an interval of interest, or having every base rendered on a
// background coloured based on its specific colour scheme. We provide many of these different
// options as 'format' methods of [`Seq`].
//
// However these helper methods themselves are constructed by building up ANSI formatted strings
// using the render module below. So if you want to implement your own functions that take a
// sequence / base and format in some specific way, this example might be for you

use seqlib::render::*;

fn main() {
    // Step 1: make any string
    let seq = String::from("ACTGCA");

    // Step 2: build up a styler with all the settings
    let styler = SeqStyler::new()
        .fixed_width()
        .colour_background(44)
        .colour_foreground(20)
        .bold();

    // Step 3: Paint sequence with the style
    let ansi_string = styler.paint(&seq);

    // Step 4: Print to any terminal that supports 8-bit colour (almost all modern terminals)
    println!("\nExample of basic styling");
    println!("{ansi_string}");

    // The real power comes when we make different stylers for different parts of our sequence
    let styler_special_base = SeqStyler::new().fixed_width().colour_background(9).bold();

    let styler_boring_bases = SeqStyler::new()
        .fixed_width()
        .colour_foreground(255)
        .colour_background(8);

    // Apply formatters to different parts of a string
    let formatted_seq = format!(
        "{}{}{}",
        styler_boring_bases.paint("ACTG"),
        styler_special_base.paint("T"),
        styler_boring_bases.paint("ACACAT")
    );

    println!("\nExample of mixing different stylers");
    println!("{formatted_seq}");

    // Can also use predefined styler themes
    println!("\nExample of predefined style (SeqStyler::HIGHLIGHT)");
    println!("{}", SeqStyler::HIGHLIGHT.paint(&seq))
}
