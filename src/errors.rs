use crate::base::Alphabet;
use thiserror::Error;

#[derive(Error, Debug)]
pub enum SeqError {
    #[error(
        "Invalid sequence for {alphabet}: found unsupported symbol(s): {invalid}. \
         Allowed symbols are standard bases plus IUPAC ambiguity codes for {alphabet}."
    )]
    InvalidCharacters { alphabet: Alphabet, invalid: String },

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
}
