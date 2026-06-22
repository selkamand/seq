use crate::{
    base::Alphabet,
    coord::{Pos, Region},
};

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum BaseError {
    #[error("invalid {alphabet} base: '{invalid}'")]
    InvalidCharacter { alphabet: Alphabet, invalid: char },

    #[error("invalid {alphabet} byte value: {invalid}; expected an ASCII nucleotide symbol")]
    InvalidByte { alphabet: Alphabet, invalid: u8 },

    #[error("base '{base}' has ambiguous chemical class")]
    AmbiguousChemicalClass { base: char },
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum CoordError {
    #[error("position must be 1-based so 0 is not a valid position")]
    PositionIsZero,

    #[error("position {value} cannot be represented on this platform; max allowed position: {max}")]
    PositionOverflowU64 { value: u64, max: Pos },

    #[error("position {value} cannot be represented on this platform; max allowed position: {max}")]
    PositionOverflowU32 { value: u32, max: Pos },

    #[error("position underflow: {lhs} - {rhs} would be < 1")]
    PositionUnderflow { lhs: Pos, rhs: usize },

    #[error("position overflow: {lhs} + {rhs} would exceed {max}")]
    PositionOverflowAdd { lhs: Pos, rhs: usize, max: Pos },

    #[error("end position of range cannot be less than start position")]
    RangeEndTooSmall { start: Pos, end: Pos },
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum SequenceError {
    #[error(transparent)]
    InvalidBase(#[from] BaseError),

    #[error(
        "invalid subsequence coordinates: requested range [{start}, {end}) on a sequence of length {len}"
    )]
    InvalidSlice {
        start: usize,
        end: usize,
        len: usize,
    },

    #[error("cannot mutate region: {region}; it spans beyond sequence length {seqlength}")]
    FailedMutateInvalidRegion { region: Region, seqlength: usize },

    #[error("cannot determine palindrome for a sequence containing ambiguous bases")]
    AmbiguousPalindrome,

    #[error(
        "cannot convert degenerate sequence to concrete sequence: ambiguous base '{base}' at position {position}"
    )]
    CannotConvertDegenerateSequence { position: Pos, base: char },
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum MutationError {
    #[error(transparent)]
    Coord(#[from] CoordError),

    #[error("mutation [{id}] does not have a sequence context")]
    MissingContext { id: String },

    #[error("cannot pyrimidine-center mutation because the reference sequence has no middle base")]
    MissingMiddleBase,

    #[error(
        "cannot pyrimidine-center mutation because middle base '{base}' has ambiguous chemical class"
    )]
    AmbiguousMiddleBase { base: char },
}

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum Error {
    #[error(transparent)]
    Base(#[from] BaseError),

    #[error(transparent)]
    Coord(#[from] CoordError),

    #[error(transparent)]
    Sequence(#[from] SequenceError),

    #[error(transparent)]
    Mutation(#[from] MutationError),
}

pub type Result<T> = std::result::Result<T, Error>;
