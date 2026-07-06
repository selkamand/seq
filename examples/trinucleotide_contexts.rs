use seqlib::coords::{Interval, Pos, Region, Strand};
use seqlib::mutations::{MutationWithContext, SmallMutation};
use seqlib::sequences::SourcedSeq;
use seqlib::{dna, pos};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let mutation = SmallMutation::new(
        "chr1".to_owned(),
        pos!(2000),
        dna!("A"),
        dna!("C"),
        Some(Strand::Positive),
        false,
        false,
    );

    let interval = Interval::new(pos!(5), pos!(10))?;

    let context = SourcedSeq::new(
        dna!("ACTGATCGAACGAGCATGCTACGGGGCCGATCGATTATCGATCAGTCA"),
        Region::new("Chr1", interval),
        Some(Strand::Positive),
    );

    let mutation_with_context = MutationWithContext::new(mutation, context, pos!(10))?;

    eprintln!("{mutation_with_context}");

    eprintln!(
        "-----Full Sequence Comparison----\n{}",
        mutation_with_context.to_difference_string()
    );

    Ok(())
}
