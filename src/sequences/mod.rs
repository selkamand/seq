mod sequence;
mod sourced;
pub use sequence::{
    BaseSliceExt, DnaSeq, IupacDnaSeq, IupacRnaSeq, RnaSeq, Seq, validate_dna_literal,
};

pub use sourced::SourcedSeq;
