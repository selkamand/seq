use seqlib::{
    coords::{Interval, Pos},
    sequences::IupacDnaSeq,
};

fn main() {
    let seq = IupacDnaSeq::new("ACTGATTTT").unwrap();

    // Highlight the 4th element in sequence vector (5th base in sequence since rust vectors are
    // 0-based)
    println!("{}", seq.format_with_highlight_index(Some(4)));

    // Highlight the 4th base by specifying its 1-based position
    let position = Pos::new(4usize).unwrap();
    println!("{}", seq.format_with_highlight_pos(Some(position)));

    // Highlight the 2nd-4th base (1-based; both-end inclusive)
    let interval = Interval::new(Pos::new(2usize).unwrap(), Pos::new(4usize).unwrap()).unwrap();
    println!("{}", seq.format_with_highlight_interval(Some(&interval)));

    // Highlight the 5th-100th base (1-based; both-end inclusive). Since seq is shorter than range,
    // will annotate with '>'
    let interval2 = Interval::new(Pos::new(5usize).unwrap(), Pos::new(100usize).unwrap()).unwrap();
    println!("{}", seq.format_with_highlight_interval(Some(&interval2)));
}
