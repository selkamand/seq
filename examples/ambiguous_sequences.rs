use seqlib::sequences::IupacDnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let clean = IupacDnaSeq::new("ACGT")?;
    let ambiguous = IupacDnaSeq::new("ACNT")?;

    assert!(clean.all_unambiguous());
    assert!(ambiguous.any_ambiguous());

    println!("Clean: {clean} (unambiguous)");
    println!("Ambiguous: {ambiguous} (contains N)");

    Ok(())
}
