use seqlib::base::Base;
use seqlib::sequence::IupacDnaSeq;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let seq = IupacDnaSeq::new("ACGT")?;

    for base in seq.as_slice() {
        println!("{}", base.to_char());
    }

    Ok(())
}
