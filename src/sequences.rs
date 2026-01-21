use crate::base::{self, Alphabet, Base};
use crate::errors::SeqError;
use core::fmt;
use std::collections::BTreeSet;

/// A biological sequence with an associated alphabet (DNA/RNA).
///
/// `Seq` represents a validated nucleotide sequence (DNA or RNA).
/// The alphabet determines which characters are considered valid and
/// influences downstream operations (e.g. reverse complementing).
///
/// This struct is not appropriate for larger-than-memory sequences.
/// It also completely ignores softmasks (for now)
#[derive(Debug, Clone)]
pub struct Seq<B: Base> {
    seq: Vec<B>, // A vector of objects with the Base Trait
    alphabet: Alphabet,
}

impl<B: Base> fmt::Display for Seq<B> {
    /// Formats the sequence as its string representation.
    ///
    /// This prints only the underlying sequence characters, without
    /// additional metadata.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use core::fmt::Write;
        for b in &self.seq {
            f.write_char(b.to_char())?;
        }
        Ok(())
    }
}

// Other functions we can run on Seq
impl<B: Base> Seq<B> {
    pub fn to_string_upper(&self) -> String {
        self.seq.iter().map(|b| b.to_char()).collect()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }
    pub fn new(sequence: &str, alphabet: Alphabet) -> Result<Seq<B>, SeqError> {
        let mut seq = Vec::with_capacity(sequence.len());

        for &byte in sequence.as_bytes() {
            // Parse one base at a time (case-insensitive handled inside try_from_ascii)
            let base = B::try_from_ascii(byte)?;
            seq.push(base);
        }

        Ok(Seq { seq, alphabet })
    }
}
// impl Seq {
//     /// Returns the length of the sequence in characters.
//     ///
//     /// This corresponds to the number of bases in the sequence.
//     pub fn len(&self) -> usize {
//         self.seq.len()
//     }
//
//     /// Returns `true` if the sequence is empty.
//     pub fn is_empty(&self) -> bool {
//         self.seq.is_empty()
//     }
//
//     /// Constructs a new `Seq` from a string slice and an alphabet.
//     ///
//     /// The sequence is validated against the provided alphabet. If any
//     /// invalid characters are found, construction fails with a
//     /// [`SeqError::InvalidCharacters`] error containing the unique set
//     /// of offending characters.
//     ///
//     /// # Errors
//     ///
//     /// Returns `SeqError::InvalidCharacters` if the sequence contains
//     /// characters not permitted by the given alphabet.
//     pub fn new(seq: &str, alphabet: Alphabet) -> Result<Self, SeqError> {
//         if seq.chars().any(|c| !alphabet.is_valid_char(c)) {
//             let invalid: BTreeSet<char> = seq
//                 .chars()
//                 .filter(|&c| !alphabet.is_valid_char(c))
//                 .collect();
//
//             if !invalid.is_empty() {
//                 let invalid = invalid.into_iter().collect::<String>();
//                 return Err(SeqError::InvalidCharacters { alphabet, invalid });
//             }
//         }
//
//         Ok(Self {
//             seq: seq.to_owned(),
//             alphabet,
//         })
//     }
//
//     /// Returns the length of the sequence in characters.
//     ///
//     /// This is an alias for [`Seq::len`].
//     pub fn length(&self) -> usize {
//         self.seq.len()
//     }
//
//     /// Returns a reference to the sequence’s alphabet.
//     pub fn alphabet(&self) -> &Alphabet {
//         &self.alphabet
//     }
// }
