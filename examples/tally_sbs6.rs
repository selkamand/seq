use seqlib::{
    base::DnaBase, coord::Pos, dna, mutation::DnaSmallMutation, pos,
    specialmuts::DnaSingleBaseSubstitution,
};

#[derive(Debug, Default)]
struct Sbs6Tally {
    c_a: u64,
    c_g: u64,
    c_t: u64,
    t_a: u64,
    t_c: u64,
    t_g: u64,
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let mutations = [
        DnaSmallMutation::new_minimal("Chr1".to_owned(), pos!(101), dna!("C"), dna!("A")),
        DnaSmallMutation::new_minimal("Chr1".to_owned(), pos!(103), dna!("C"), dna!("T")),
        DnaSmallMutation::new_minimal("Chr1".to_owned(), pos!(104), dna!("T"), dna!("A")),
        DnaSmallMutation::new_minimal("Chr1".to_owned(), pos!(105), dna!("T"), dna!("C")),
        DnaSmallMutation::new_minimal("Chr1".to_owned(), pos!(106), dna!("T"), dna!("G")),
        DnaSmallMutation::new_minimal("Chr1".to_owned(), pos!(107), dna!("G"), dna!("A")),
        DnaSmallMutation::new_minimal("Chr1".to_owned(), pos!(108), dna!("C"), dna!("C")),
        DnaSmallMutation::new_minimal("Chr1".to_owned(), pos!(109), dna!("C"), dna!("CA")),
    ];

    let mut tally = Sbs6Tally::default();

    for mutation in &mutations {
        // Parse single base mutation from general small mutation format
        let sbs = match DnaSingleBaseSubstitution::try_from(mutation) {
            Ok(sbs) => sbs,
            Err(err) => {
                eprintln!(
                    "Skipping mutation: {} as it was not a valid single base substition: {err}",
                    mutation.chrom_pos_ref_alt(),
                );
                continue;
            }
        };

        let sbs_pyrimidine_centered = match sbs.try_pyrimidine_center() {
            Ok(s) => s,
            Err(err) => {
                eprintln!(
                    "Skipping mutation: {} as it was not a valid single base substition: {err}",
                    mutation.chrom_pos_ref_alt()
                );
                continue;
            }
        };

        match (
            *sbs_pyrimidine_centered.reference(),
            *sbs_pyrimidine_centered.alternative(),
        ) {
            (DnaBase::C, DnaBase::A) => tally.c_a += 1,
            (DnaBase::C, DnaBase::G) => tally.c_g += 1,
            (DnaBase::C, DnaBase::T) => tally.c_t += 1,
            (DnaBase::T, DnaBase::A) => tally.t_a += 1,
            (DnaBase::T, DnaBase::C) => tally.t_c += 1,
            (DnaBase::T, DnaBase::G) => tally.t_g += 1,
            _ => unreachable!("Pyrimidine centering does not allow non-C/T reference bases"),
        }
    }

    println!("{tally:#?}");

    Ok(())
}
