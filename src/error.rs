use crate::{
    base::Alphabet,
    coord::{Pos, Region},
    mutation::SmallMutation,
};

#[derive(thiserror::Error, Debug, PartialEq, Eq)]
pub enum Error {
    #[error(
        "Invalid {alphabet} base: '{invalid}'. \
         Allowed symbols are standard bases plus IUPAC ambiguity codes for {alphabet}."
    )]
    InvalidCharacter { alphabet: Alphabet, invalid: char },

    #[error(
        "Invalid {alphabet} byte value: {invalid}. \
         Expected an ASCII letter representing a nucleotide (e.g. A,C,G,T/U,N)."
    )]
    InvalidByte { alphabet: Alphabet, invalid: u8 },

    #[error(
        "Invalid subsequence coordinates: requested range [{start}, {end}) on a sequence of length {len}. \
         Indices are 0-based and the end position is exclusive (like Rust slicing)."
    )]
    InvalidSlice {
        start: usize,
        end: usize,
        len: usize,
    },

    #[error("Position must be 1-based so 0 is not a valid position")]
    PositionIsZero,

    #[error("position {value} cannot be represented on this platform. Max allowed position: {max}")]
    PositionOverflowU64 { value: u64, max: Pos },

    #[error("position {value} cannot be represented on this platform. Max allowed position: {max}")]
    PositionOverflowU32 { value: u32, max: Pos },

    #[error("position underflow: {lhs} - {rhs} would be < 1")]
    PositionUnderflow { lhs: Pos, rhs: usize },

    #[error("position overflow: {lhs} + {rhs} would exceed {max}")]
    PositionOverflowAdd { lhs: Pos, rhs: usize, max: Pos },

    #[error("End position of range: [end] cannot be less than start position [start]")]
    RangeEndTooSmall { start: Pos, end: Pos },

    #[error("Mutation [{id}] does not have a sequence context")]
    MutationMissingContext { id: String },

    #[error("Can NOT mutate region: {region}. Spans beyond sequence ({seqlength} bases long")]
    FailedMutateInvalidRegion { region: Region, seqlength: usize },

    #[error("Can NOT determine palindrome: {0}")]
    PalindromeError(PalindromeErrorReason),

    #[error(
        "Cannot convert degenerate sequence to concrete sequence: ambiguous base '{base}' at position {position}"
    )]
    CannotConvertDegenerateSequence { position: Pos, base: char },

    #[error(
        "Cannot pyrmidine-center the mutation. Either sequence was either length or middlebase was a degenerate Iupac character (like N)"
    )]
    InvalidMiddlebaseCannotPyrimidineCenter,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum PalindromeErrorReason {
    AmbiguousBases,
}
impl std::fmt::Display for PalindromeErrorReason {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            PalindromeErrorReason::AmbiguousBases => write!(f, "AmbiguousBases"),
        }
    }
}

pub type Result<T> = std::result::Result<T, Error>;
