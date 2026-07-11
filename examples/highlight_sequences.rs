use seqlib::{
    coords::{Interval, Pos},
    pos,
    sequences::IupacDnaSeq,
};

fn main() {
    let seq = IupacDnaSeq::new("ACTGATTTT").unwrap();

    // Highlight the 4th element in sequence vector (5th base in sequence since rust vectors are
    // 0-based)
    println!("Highlight the 4th base in the sequence");
    let position = pos!(4);
    println!("{}", seq.format_with_highlight_pos(Some(position)));

    // Highlight the 2nd-4th base (1-based; both-end inclusive)
    println!("\nHighlight the interval: 2nd-4th (1-based both-end inclusive) with square brakcets");
    let interval = Interval::new(pos!(2), pos!(4)).unwrap();
    println!("{}", seq.format_with_highlight_interval(Some(&interval)));

    println!("\nHighlight the interval: 2nd-4th (1-based both-end inclusive) using ansi colour");
    println!("{}", seq.format_with_coloured_interval(&interval));

    // Highlight the 5th-100th base (1-based; both-end inclusive). Since seq is shorter than range,
    // will annotate with '>'
    println!("\nHighlight the 5th-100th base with text:");
    let interval2 = Interval::new(Pos::new(5usize).unwrap(), Pos::new(100usize).unwrap()).unwrap();
    println!("{}", seq.format_with_highlight_interval(Some(&interval2)));

    println!("\nHighlight the 5th-100th base with ansi colours:");
    println!("{}", seq.format_with_coloured_interval(&interval2));

    // Pretty print with background colours
    println!("\nColour bases using ansi:");
    println!("{}", seq.format_with_colour());
}
