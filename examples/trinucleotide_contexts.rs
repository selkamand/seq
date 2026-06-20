use seqlib::coord::{Pos, Strand};
use seqlib::{
    context::{self, Orientation},
    error::Result,
    mutation::{MutationWithContext, SmallMutation},
};
use seqlib::{dna, pos};

fn main() -> Result<()> {
    let mutation = SmallMutation::new(
        "chr1".to_owned(),
        pos!(2000),
        dna!("A"),
        dna!("C"),
        Some(Strand::Positive),
        false,
        false,
    );

    let context = context::ContextWindow::new(
        dna!("ACTGATCGATCGAGCATGCTACGGGGCCGATCGATTATCGATCAGTCA"),
        pos!(10),
        pos!(5),
        Orientation::Forward,
    );

    let mutation_with_context = MutationWithContext::new(mutation, Some(context));

    println!("{mutation_with_context}");

    Ok(())
}
