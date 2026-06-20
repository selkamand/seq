//! Special mutation types
//!
//! Small mutations can all be represented by the [`SmallMutation`] types, but certain downstream
//! operations are much simpler to perform on minimal, highly specialised structures.
//! For example pyrimidine centering and subsequent classification of single base substitution types
//! for mutational signature analysis is simpler if instead of having reference and alternative 'sequences' you have just the individual
//! bases

use crate::{
    base::Base,
    mutation::{SmallMutation, SmallMutationType},
};

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum Error<B: Base> {
    #[error(". Mutation was probably not a SNV {mutation}")]
    TryFromSmallMutationErrorWrongClass {
        mutation: SmallMutation<B>,
        class: SmallMutationType,
    },

    #[error("Operation does not support ambiguous base: {base}")]
    InvalidBaseAmbiguous { base: B },
}

/// A simple representation of a single base substition
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub struct SingleBaseSubstition<B: Base> {
    reference: B,
    alternative: B,
}

impl<B: Base> SingleBaseSubstition<B> {
    /// Reverse complement a single base substitution
    /// Complements the reference and alternative bases. Note we don't need to reverse
    /// because in a single base subtition there is only one nucleotide in reference or alt
    pub fn reverse_complement(&self) -> Self {
        Self {
            reference: self.alternative.complement(),
            alternative: self.reference.complement(),
        }
    }

    ///  Ensure reference base is a pyrimidine
    ///
    /// This is accomplished by reverse complementing the mutation
    /// if reference is a purine.
    ///
    /// # Errors
    /// If the chemical class of the mutation is ambiguous (e.g. because of the use of Iupac Bases
    /// that support degenerate bases) we return an [`Error::InvalidBaseAmbiguous`] error
    ///
    pub fn pyrimidine_center(&self) -> Result<Self, Error<B>> {
        match self.reference.chemical_class() {
            crate::base::ChemClass::Purine => Ok(self.reverse_complement()),
            crate::base::ChemClass::Pyrimidine => Ok(*self),
            crate::base::ChemClass::Ambiguous => Err(Error::InvalidBaseAmbiguous {
                base: self.reference,
            }),
        }
    }
}

// Convert SmallMutation to Single Base Substitution
impl<B: Base> TryFrom<SmallMutation<B>> for SingleBaseSubstition<B> {
    type Error = Error<B>;

    fn try_from(value: SmallMutation<B>) -> Result<Self, Self::Error> {
        let class = value.class();

        match class {
            SmallMutationType::SNV => {
                let reference = value.reference().as_slice()[0];
                let alternative = value.alternative().as_slice()[0];

                Ok(Self {
                    reference,
                    alternative,
                })
            }

            other => Err(Error::TryFromSmallMutationErrorWrongClass {
                mutation: value,
                class: other,
            }),
        }
    }
}
