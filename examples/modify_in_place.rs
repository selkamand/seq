/// Sequence transformations like complement, reverse_complement, etc.
/// all create a new owned copy of the original with the transformations applied.
/// The original Seq object is unmodified.
///
/// This is great for usability but sometimes you may want to modify the original object to
/// avoid creating a whole new one in memory.
///
/// For this reason, we provide '_in_place' versions of all transformations that modify the
/// original struct (and return nothing).

fn main() -> Result<(), seqlib::errors::SeqError> {
    // Create mutable sequence
    let mut seq = seqlib::sequences::DnaSeq::new("AGACT").unwrap();
    eprintln!("original: {seq}");

    // Complement(modify in place)
    seq.complement_in_place();
    eprintln!("complement: {seq}");

    // Complement back
    seq.complement_in_place();
    eprintln!("complimented again: {seq}");

    Ok(())
}
