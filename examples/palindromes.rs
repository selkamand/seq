use seqlib::sequences::{DnaSeq, IupacDnaSeq};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let palindrome = DnaSeq::new("GAATTC")?;
    let ambiguous = IupacDnaSeq::new("NNNNNN")?;

    assert!(palindrome.is_palindromic());
    assert!(ambiguous.is_palindromic_checked().is_err());

    println!("{palindrome} is palindromic");
    println!("{ambiguous} is NOT palindromic (ambiguous)");

    Ok(())
}
