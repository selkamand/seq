use core::fmt;

use crate::{
    base::{Base, ChemClass},
    sequences::DnaSeq,
};

// Small Variant Class
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmallMutation {
    chromosome: String,
    position: u64, // 1-based position (start)
    reference: DnaSeq,
    alternative: DnaSeq,
    multiallelic: bool,
    context: Option<DnaSeq>,
}

// Implement the `fmt::Display` trait for `Point`.
impl fmt::Display for SmallMutation {
    // This trait requires the `fmt` function with this exact signature.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use the write! macro to format the output.
        write!(
            f,
            "{}:{} {}>{} (delta: {}; class: {}; multiallelic:{})",
            self.chromosome,
            self.position,
            self.reference,
            self.alternative,
            self.delta(),
            self.class(),
            self.multiallelic
        )
    }
}

impl SmallMutation {
    pub fn new(
        chromosome: String,
        position: u64,
        reference: DnaSeq,
        alternative: DnaSeq,
        multiallelic: bool,
        context: Option<DnaSeq>,
    ) -> Self {
        Self {
            chromosome,
            position,
            reference,
            alternative,
            multiallelic,
            context,
        }
    }

    pub fn add_context(&mut self, seq: DnaSeq) {
        self.context = Some(seq)
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmallMutationType {
    SNV,
    DOUBLET,
    MNV,
    INSERTION,
    DELETION,
}

impl fmt::Display for SmallMutationType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            SmallMutationType::SNV => "SNV",
            SmallMutationType::DOUBLET => "DOUBLET",
            SmallMutationType::MNV => "MNV",
            SmallMutationType::INSERTION => "INSERTION",
            SmallMutationType::DELETION => "DELETION",
        };
        write!(f, "{s}")
    }
}

impl SmallMutationType {
    pub fn from_lengths(reflen: usize, altlen: usize) -> Self {
        match altlen.cmp(&reflen) {
            std::cmp::Ordering::Greater => Self::INSERTION,
            std::cmp::Ordering::Less => Self::DELETION,
            std::cmp::Ordering::Equal => match reflen {
                0 => Self::MNV, // OR HANDLE AS ERROR ELSEWHERE (0-LENGTH ALLELES ARE INVALID FOR VCF)
                1 => Self::SNV,
                2 => Self::DOUBLET,
                _ => Self::MNV,
            },
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TiTv {
    Transition,
    Transversion,
}

impl TiTv {
    pub fn from_chemical_class(reference: ChemClass, alternative: ChemClass) -> Option<TiTv> {
        match (reference, alternative) {
            (ChemClass::Purine, ChemClass::Purine) => Some(TiTv::Transition),
            (ChemClass::Pyrimidine, ChemClass::Pyrimidine) => Some(TiTv::Transition),
            (ChemClass::Purine, ChemClass::Pyrimidine) => Some(TiTv::Transversion),
            (ChemClass::Pyrimidine, ChemClass::Purine) => Some(TiTv::Transversion),
            _ => None, // If either chemical class is ambiguous, return None
        }
    }
}

impl SmallMutation {
    pub fn reflen(&self) -> usize {
        self.reference.len()
    }
    pub fn altlen(&self) -> usize {
        self.alternative.len()
    }

    pub fn delta(&self) -> i64 {
        self.altlen() as i64 - self.reflen() as i64
    }

    pub fn class(&self) -> SmallMutationType {
        SmallMutationType::from_lengths(self.reflen(), self.altlen())
    }

    pub fn titv(&self) -> Option<TiTv> {
        if self.class() != SmallMutationType::SNV {
            return None;
        }

        let r = self.reference.as_slice().first()?;
        let a = self.alternative.as_slice().first()?;

        TiTv::from_chemical_class(r.chemical_class(), a.chemical_class())
    }
}
