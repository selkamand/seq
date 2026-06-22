//! Special mutation types
//!
//! Small mutations (Single base substitions, Doublets, MNVs, INDELs)can all be represented by the [`SmallMutation`] types, but certain downstream
//! operations are much simpler to perform on minimal, highly specialised structures.
//! For example pyrimidine centering and subsequent classification of single base substitution types
//! for mutational signature analysis is simpler if instead of having reference and alternative 'sequences'
//! (vectors of bases).
//! you force ref and alt to be individual bases
use crate::{
    base::{Base, DegenerateBase, DnaBase, IupacDnaBase, IupacRnaBase, RnaBase},
    mutation::{SmallMutation, SmallMutationType},
};

/// Error raised when the [`SmallMutation`] types cannot be narrowed to a concrete SBS.
#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum SingleBaseSubstitutionError {
    #[error("mutation class {class} is not a single-nucleotide variant")]
    WrongClass { class: SmallMutationType },

    #[error(
        "reference base '{base}' is ambiguous so can NOT be represented as a single concrete base"
    )]
    AmbiguousReference { base: char },

    #[error(
        "alternative base '{base}' is ambiguous so can NOT be represented as a single, concrete base"
    )]
    AmbiguousAlternative { base: char },
}

/// A simple representation of a single base substition
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct SingleBaseSubstitution<B: Base> {
    pub reference: B,
    pub alternative: B,
}

// Add pub types for concrete base sequences
pub type DnaSingleBaseSubstitution = SingleBaseSubstitution<DnaBase>;
pub type RnaSingleBaseSubstitution = SingleBaseSubstitution<RnaBase>;

impl<B: Base> SingleBaseSubstitution<B> {
    /// Constructor
    pub fn new(reference: B, alternative: B) -> Self {
        Self {
            reference,
            alternative,
        }
    }

    /// Reverse complement a single base substitution
    /// Complements the reference and alternative bases. Note we don't need to reverse
    /// because in a single base subtition there is only one nucleotide in reference or alt
    pub fn reverse_complement(&self) -> Self {
        Self {
            reference: self.alternative.complement(),
            alternative: self.reference.complement(),
        }
    }

    /// Getter for reference base
    pub fn reference(&self) -> &B {
        &self.reference
    }

    /// Getter for alternative base
    pub fn alternative(&self) -> &B {
        &self.alternative
    }
    /// Ensure the reference base is a pyrimidine.
    ///
    /// If the reference base is a purine, the substitution is reverse-complemented.
    ///
    /// # Errors
    ///
    /// Returns [`SingleBaseSubstitutionError::AmbiguousReference`] if the reference base
    /// does not have a certain purine/pyrimidine class.
    pub fn pyrimidine_center(&self) -> Result<Self, SingleBaseSubstitutionError> {
        match self.reference.chemical_class() {
            crate::base::ChemClass::Purine => Ok(self.reverse_complement()),
            crate::base::ChemClass::Pyrimidine => Ok(*self),
            crate::base::ChemClass::Ambiguous => {
                Err(SingleBaseSubstitutionError::AmbiguousReference {
                    base: self.reference.to_char(),
                })
            }
        }
    }
}

// Convert SmallMutation to Single Base Substitution when base type matches
impl<B: Base> TryFrom<&SmallMutation<B>> for SingleBaseSubstitution<B> {
    type Error = SingleBaseSubstitutionError;

    fn try_from(value: &SmallMutation<B>) -> Result<Self, Self::Error> {
        ensure_snv(value)?;
        Ok(Self::new(
            value.reference().as_slice()[0],
            value.alternative().as_slice()[0],
        ))
    }
}

impl TryFrom<&SmallMutation<IupacDnaBase>> for DnaSingleBaseSubstitution {
    type Error = SingleBaseSubstitutionError;

    fn try_from(value: &SmallMutation<IupacDnaBase>) -> Result<Self, Self::Error> {
        ensure_snv(value)?;

        let reference = value.reference().as_slice()[0]
            .try_to_concrete()
            .ok_or_else(|| SingleBaseSubstitutionError::AmbiguousReference {
                base: value.reference().as_slice()[0].to_char(),
            })?;

        let alternative = value.alternative().as_slice()[0]
            .try_to_concrete()
            .ok_or_else(|| SingleBaseSubstitutionError::AmbiguousAlternative {
                base: value.alternative().as_slice()[0].to_char(),
            })?;

        Ok(Self::new(reference, alternative))
    }
}

impl TryFrom<&SmallMutation<IupacRnaBase>> for RnaSingleBaseSubstitution {
    type Error = SingleBaseSubstitutionError;

    fn try_from(value: &SmallMutation<IupacRnaBase>) -> Result<Self, Self::Error> {
        ensure_snv(value)?;

        let reference = value.reference().as_slice()[0]
            .try_to_concrete()
            .ok_or_else(|| SingleBaseSubstitutionError::AmbiguousReference {
                base: value.reference().as_slice()[0].to_char(),
            })?;
        let alternative = value.alternative().as_slice()[0]
            .try_to_concrete()
            .ok_or_else(|| SingleBaseSubstitutionError::AmbiguousAlternative {
                base: value.alternative().as_slice()[0].to_char(),
            })?;

        Ok(Self::new(reference, alternative))
    }
}

fn ensure_snv<B: Base>(mutation: &SmallMutation<B>) -> Result<(), SingleBaseSubstitutionError> {
    let class = mutation.class();
    if class == SmallMutationType::SNV {
        Ok(())
    } else {
        Err(SingleBaseSubstitutionError::WrongClass { class })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        coord::Pos,
        mutation::IupacDnaSmallMutation,
        sequence::{DnaSeq, IupacDnaSeq},
    };

    fn iupac_dna_mut(reference: &str, alternative: &str) -> IupacDnaSmallMutation {
        IupacDnaSmallMutation::new(
            "chr1".to_string(),
            Pos::new(1).unwrap(),
            IupacDnaSeq::new(reference).unwrap(),
            IupacDnaSeq::new(alternative).unwrap(),
            None,
            false,
            true,
        )
    }

    #[test]
    fn converts_iupac_snv_to_concrete_dna_sbs() {
        let mutation = iupac_dna_mut("C", "A");
        let sbs = DnaSingleBaseSubstitution::try_from(&mutation).unwrap();

        assert_eq!(sbs.reference(), &DnaBase::C);
        assert_eq!(sbs.alternative(), &DnaBase::A);
    }

    #[test]
    fn rejects_wrong_class() {
        let mutation = iupac_dna_mut("CC", "AA");
        let err = DnaSingleBaseSubstitution::try_from(&mutation).unwrap_err();

        assert_eq!(
            err,
            SingleBaseSubstitutionError::WrongClass {
                class: SmallMutationType::DOUBLET
            }
        );
    }

    #[test]
    fn rejects_ambiguous_reference_and_alternative_separately() {
        let ref_err = DnaSingleBaseSubstitution::try_from(&iupac_dna_mut("N", "A")).unwrap_err();
        assert_eq!(
            ref_err,
            SingleBaseSubstitutionError::AmbiguousReference { base: 'N' }
        );

        let alt_err = DnaSingleBaseSubstitution::try_from(&iupac_dna_mut("C", "N")).unwrap_err();
        assert_eq!(
            alt_err,
            SingleBaseSubstitutionError::AmbiguousAlternative { base: 'N' }
        );
    }

    #[test]
    fn pyrimidine_center_reverse_complements_purine_reference() {
        let sbs = DnaSingleBaseSubstitution::new(DnaBase::G, DnaBase::A);
        let centered = sbs.pyrimidine_center().unwrap();

        assert_eq!(centered.reference(), &DnaBase::T);
        assert_eq!(centered.alternative(), &DnaBase::C);
    }

    #[test]
    fn converts_concrete_dna_small_mutation() {
        let mutation = crate::mutation::DnaSmallMutation::new(
            "chr1".to_string(),
            Pos::new(1).unwrap(),
            DnaSeq::new("T").unwrap(),
            DnaSeq::new("G").unwrap(),
            None,
            false,
            true,
        );

        let sbs = DnaSingleBaseSubstitution::try_from(&mutation).unwrap();
        assert_eq!(sbs.reference(), &DnaBase::T);
        assert_eq!(sbs.alternative(), &DnaBase::G);
    }
}
