use seqlib::sequence::IupacDnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let palindrome = IupacDnaSeq::new("GAATTC")?;
    let ambiguous = IupacDnaSeq::new("NNNNNN")?;

    assert!(palindrome.is_palindromic());
    assert!(!ambiguous.is_palindromic());

    println!("{palindrome} is palindromic");
    println!("{ambiguous} is NOT palindromic (ambiguous)");

    Ok(())
}
