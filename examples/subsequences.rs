use seqlib::{
    coord::{Pos, Region},
    pos,
    sequence::{BaseSliceExt, DnaSeq},
};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    // Define a sequence
    let seq = DnaSeq::new("ACGTAC")?;
    println!("{seq} <- Sequence (original)");

    // Define a Region (1 based start & end, both-end inclusive)
    let start = pos!(2);
    let end = pos!(4);
    let region = Region::new(start, end)?;

    // Highlight where this region is on our sequence
    println!(
        "{} <- Sequence (annotated by region {} [{}bp])",
        seq.format_with_highlight_region(Some(&region)),
        region,
        region.len()
    );

    // Grab the subsequence (owned copy)
    let subseq = seq.subseq(&region)?;

    // Print out the slice with
    println!("{subseq} <- sub-sequence");

    // If you just want to borrow a slice, use the subseq_slice method
    let subseq_slice = seq.subseq_slice(&region)?;
    println!("{} <- sub-sequence", subseq_slice.to_string_upper());

    Ok(())
}

