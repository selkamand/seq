use seqlib::coords::{Interval1, Pos1, Region, Strand};
use seqlib::mutations::{MutationWithContext, SmallMutation};
use seqlib::sequences::{BaseSliceExt, SourcedSeq};
use seqlib::{dna, pos1};
use std::error::Error;

fn main() -> Result<(), Box<dyn Error>> {
    let mutation = SmallMutation::new(
        "Chr1".to_owned(),
        pos1!(2000),
        dna!("A"),
        dna!("C"),
        Some(Strand::Positive),
    );

    let interval = Interval1::new(pos1!(2000), pos1!(2000))?;

    let context = SourcedSeq::new(
        dna!("ACTGATCGAACGAGCATGCTACGGGGCCGATCGATTATCGATCAGTCA"),
        Region::new("Chr1", interval),
        Some(Strand::Positive),
    );

    let mutation_with_context = MutationWithContext::new(mutation, context)?;

    eprintln!("{mutation_with_context}");

    eprintln!(
        "-----Full Sequence Comparison----\n{}",
        mutation_with_context.to_difference_string()
    );

    let tnc = mutation_with_context.kmer_centered_on_anchor(3);
    eprintln!(
        "\n\nTNC: {}",
        tnc.map(|x| x.to_string_upper())
            .unwrap_or("Could not be found".to_string())
    );

    Ok(())
}
