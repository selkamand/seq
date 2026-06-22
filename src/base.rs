// Data structures for representing and transforming individual nucleotides.
// We currently implement DnaBase and RnaBase, which both implement the Base trait.
// The Base trait defines the small set of operations any "nucleotide alphabet" must support.

use core::fmt;

use crate::error::BaseError;

/// A minimal interface for “a single nucleotide symbol” (DNA, RNA, or a future alphabet).
///
/// Think of this trait as: *"What can I do with one base?"*
///
/// If you implement `Base` for an enum, you get:
/// - a **complement** mapping (A<->T/U, ambiguity codes too)
/// - conversion to **ASCII** bytes for printing / writing to files
/// - parsing from an ASCII byte with helpful errors
/// - a flag telling you if the base is ambiguous (e.g. R means A/G)
///
/// Why ASCII bytes (`u8`) instead of `char`?
/// - Real sequence files are bytes.
/// - DNA/RNA/IUPAC codes are ASCII characters.
/// - Using bytes is fast and avoids Unicode complexity.
///
/// The contract is:
/// - `complement` is infallible (it always returns something)
/// - `try_from_ascii` is the validation gate (it may fail)
pub trait Base: Copy + Eq + fmt::Debug + fmt::Display + Sized {
    /// Name of the alphabet
    const ALPHABET: Alphabet;

    /// Return the complement of this base.
    ///
    /// Examples (DNA):
    /// - A ↔ T
    /// - C ↔ G
    /// - R (A/G) ↔ Y (C/T)
    ///
    /// This must be infallible: every valid base must have a defined complement.
    fn complement(self) -> Self;

    /// Convert this base to an uppercase ASCII byte (e.g. `b'A'`).
    ///
    /// This is useful for:
    /// - printing
    /// - writing FASTA/FASTQ
    /// - building `String`s
    fn to_ascii(self) -> u8;

    fn to_char(self) -> char {
        self.to_ascii() as char
    }

    /// Convert this base to a lowercase ASCII byte (e.g. `b'a'`).
    ///
    /// Lowercase bases are often used for “soft-masking” (e.g. low-confidence regions),
    /// even though the biological base is the same.
    fn to_ascii_lower(self) -> u8;

    /// Check if base only represents one possible nucleotide (e.g. A/C/T/U/G).
    /// Bases like `R` can represent multiple possible nucleotides (in this case, A or G)
    fn is_unambiguous(self) -> bool {
        !self.is_ambiguous()
    }

    /// Parse a single ASCII byte into a base.
    ///
    /// This is the “gatekeeper” function: it checks whether an input symbol is allowed.
    ///
    /// Rules:
    /// - input must be ASCII (bytes 0–127)
    /// - parsing is case-insensitive (both `b'a'` and `b'A'` work)
    /// - invalid input returns a `SeqError` describing what went wrong
    fn try_from_ascii(b: u8) -> Result<Self, BaseError>;

    /// Returns `true` if this symbol can represent more than one concrete nucleotide.
    ///
    /// Examples:
    /// - `A` is not ambiguous (it always means A)
    /// - `R` is ambiguous (it means A **or** G)
    /// - `N` is ambiguous/unknown (it means any base)
    ///
    /// This is useful because some operations (like translation) usually require
    /// unambiguous sequences.
    fn is_ambiguous(self) -> bool;

    /// Classify the type of Base (Purine Vs Pyrimidine)
    ///
    /// # Error
    /// Returns [`BaseError::AmbiguousChemicalClass`] if the chemical class is unclear
    /// (happens when iupac base is ambiguous)
    fn try_chemical_class(self) -> Result<ChemClass, BaseError>;
}

/// Chemical class of a nucleotide base (purine/pyrimidine).
///
/// This enum describes whether a base:
/// - is a purine,
/// - is a pyrimidine
///
/// Examples:
/// - `A` / G` → `Purine`
/// - `T` / `C` / `U`  → `Pyrimidine`
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum ChemClass {
    /// The base is a purine (A or G, or ambiguity codes that
    /// only expand to purines, e.g. `R`).
    Purine,

    /// The base is certainly a pyrimidine (C or T/U, or ambiguity codes that
    /// only expand to pyrimidines, e.g. `Y`).
    Pyrimidine,
}

/// DNA nucleotide symbols including IUPAC ambiguity codes.
///
/// Stored as an enum so Rust can enforce that functions like `complement` handle
/// every possible symbol (no missing cases).
#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum IupacDnaBase {
    A,
    C,
    G,
    T,
    /// Unknown / any base
    N,
    /// A or G
    R,
    /// C or T
    Y,
    /// G or C
    S,
    /// A or T
    W,
    /// G or T
    K,
    /// A or C
    M,
    /// C or G or T
    B,
    /// A or G or T
    D,
    /// A or C or T
    H,
    /// A or C or G
    V,
}

/// RNA nucleotide symbols including IUPAC ambiguity codes.
#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum IupacRnaBase {
    A,
    C,
    G,
    U,
    /// Unknown / any base
    N,
    R,
    Y,
    S,
    W,
    K,
    M,
    B,
    D,
    H,
    V,
}

/// DNA nucleotide symbols not including IUPAC ambiguity codes.
#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum DnaBase {
    A,
    C,
    G,
    T,
}

/// RNA nucleotide symbols not including IUPAC ambiguity codes.
#[repr(u8)]
#[derive(Clone, Copy, Debug, Eq, PartialEq, Ord, PartialOrd)]
pub enum RnaBase {
    A,
    C,
    G,
    U,
}

impl Base for IupacDnaBase {
    const ALPHABET: Alphabet = Alphabet::DNA;

    fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::N => Self::N,
            Self::R => Self::Y,
            Self::Y => Self::R,
            Self::S => Self::S,
            Self::W => Self::W,
            Self::K => Self::M,
            Self::M => Self::K,
            Self::B => Self::V,
            Self::V => Self::B,
            Self::D => Self::H,
            Self::H => Self::D,
        }
    }

    fn try_from_ascii(b: u8) -> Result<Self, BaseError> {
        if !b.is_ascii() {
            return Err(BaseError::InvalidByte {
                alphabet: Self::ALPHABET,
                invalid: b,
            });
        }

        Self::from_ascii_const(b).ok_or(BaseError::InvalidCharacter {
            alphabet: Self::ALPHABET,
            invalid: b as char,
        })
    }

    fn to_ascii(self) -> u8 {
        match self {
            IupacDnaBase::A => b'A',
            IupacDnaBase::C => b'C',
            IupacDnaBase::G => b'G',
            IupacDnaBase::T => b'T',
            IupacDnaBase::N => b'N',
            IupacDnaBase::R => b'R',
            IupacDnaBase::Y => b'Y',
            IupacDnaBase::S => b'S',
            IupacDnaBase::W => b'W',
            IupacDnaBase::K => b'K',
            IupacDnaBase::M => b'M',
            IupacDnaBase::B => b'B',
            IupacDnaBase::D => b'D',
            IupacDnaBase::H => b'H',
            IupacDnaBase::V => b'V',
        }
    }

    fn to_ascii_lower(self) -> u8 {
        match self {
            IupacDnaBase::A => b'a',
            IupacDnaBase::C => b'c',
            IupacDnaBase::G => b'g',
            IupacDnaBase::T => b't',
            IupacDnaBase::N => b'n',
            IupacDnaBase::R => b'r',
            IupacDnaBase::Y => b'y',
            IupacDnaBase::S => b's',
            IupacDnaBase::W => b'w',
            IupacDnaBase::K => b'k',
            IupacDnaBase::M => b'm',
            IupacDnaBase::B => b'b',
            IupacDnaBase::D => b'd',
            IupacDnaBase::H => b'h',
            IupacDnaBase::V => b'v',
        }
    }

    fn is_ambiguous(self) -> bool {
        !matches!(
            self,
            IupacDnaBase::A | IupacDnaBase::C | IupacDnaBase::G | IupacDnaBase::T
        )
    }

    fn try_chemical_class(self) -> Result<ChemClass, BaseError> {
        match self {
            IupacDnaBase::A => Ok(ChemClass::Purine),
            IupacDnaBase::G => Ok(ChemClass::Purine),
            IupacDnaBase::C => Ok(ChemClass::Pyrimidine),
            IupacDnaBase::T => Ok(ChemClass::Pyrimidine),
            IupacDnaBase::R => Ok(ChemClass::Purine),
            IupacDnaBase::Y => Ok(ChemClass::Pyrimidine),
            IupacDnaBase::N => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::N.to_char(),
            }),
            IupacDnaBase::S => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::S.to_char(),
            }),
            IupacDnaBase::W => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::W.to_char(),
            }),
            IupacDnaBase::K => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::K.to_char(),
            }),
            IupacDnaBase::M => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::M.to_char(),
            }),
            IupacDnaBase::B => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::B.to_char(),
            }),
            IupacDnaBase::D => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::D.to_char(),
            }),
            IupacDnaBase::H => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::H.to_char(),
            }),
            IupacDnaBase::V => Err(BaseError::AmbiguousChemicalClass {
                base: IupacDnaBase::V.to_char(),
            }),
        }
    }
}

impl fmt::Display for IupacDnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use core::fmt::Write;
        f.write_char((*self).to_char())
    }
}
impl fmt::Display for DnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use core::fmt::Write;
        f.write_char((*self).to_char())
    }
}
impl fmt::Display for RnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use core::fmt::Write;
        f.write_char((*self).to_char())
    }
}

impl fmt::Display for IupacRnaBase {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        use core::fmt::Write;
        f.write_char((*self).to_char())
    }
}

impl Base for DnaBase {
    const ALPHABET: Alphabet = Alphabet::DNA;

    fn complement(self) -> Self {
        match self {
            Self::A => Self::T,
            Self::T => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
        }
    }

    fn try_from_ascii(b: u8) -> Result<Self, BaseError> {
        if !b.is_ascii() {
            return Err(BaseError::InvalidByte {
                alphabet: Self::ALPHABET,
                invalid: b,
            });
        }

        Self::from_ascii_const(b).ok_or(BaseError::InvalidCharacter {
            alphabet: Self::ALPHABET,
            invalid: b as char,
        })
    }

    fn to_ascii(self) -> u8 {
        match self {
            DnaBase::A => b'A',
            DnaBase::C => b'C',
            DnaBase::G => b'G',
            DnaBase::T => b'T',
        }
    }

    fn to_ascii_lower(self) -> u8 {
        match self {
            DnaBase::A => b'a',
            DnaBase::C => b'c',
            DnaBase::G => b'g',
            DnaBase::T => b't',
        }
    }

    fn is_ambiguous(self) -> bool {
        false
    }

    fn try_chemical_class(self) -> Result<ChemClass, BaseError> {
        match self {
            DnaBase::A => Ok(ChemClass::Purine),
            DnaBase::G => Ok(ChemClass::Purine),
            DnaBase::C => Ok(ChemClass::Pyrimidine),
            DnaBase::T => Ok(ChemClass::Pyrimidine),
        }
    }
}

impl IupacDnaBase {
    /// Parse a single ASCII byte into a DNA base (A,C,G,T plus IUPAC codes),
    /// case-insensitive.
    ///
    /// Returns `None` if `b` is not a valid DNA symbol.
    pub const fn from_ascii_const(b: u8) -> Option<Self> {
        match b {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'T' | b't' => Some(Self::T),

            b'N' | b'n' => Some(Self::N),
            b'R' | b'r' => Some(Self::R),
            b'Y' | b'y' => Some(Self::Y),
            b'S' | b's' => Some(Self::S),
            b'W' | b'w' => Some(Self::W),
            b'K' | b'k' => Some(Self::K),
            b'M' | b'm' => Some(Self::M),
            b'B' | b'b' => Some(Self::B),
            b'D' | b'd' => Some(Self::D),
            b'H' | b'h' => Some(Self::H),
            b'V' | b'v' => Some(Self::V),

            _ => None,
        }
    }
}

impl IupacRnaBase {
    /// Parse a single ASCII byte into an [`RnaBase`] in a `const` context.
    ///
    /// This is a const-friendly equivalent of [`Base::try_from_ascii`], intended
    /// for compile-time validation (e.g. in macros).
    ///
    /// - Accepts both uppercase and lowercase ASCII letters
    /// - Returns `None` for invalid characters
    /// - Does not panic
    pub const fn from_ascii_const(b: u8) -> Option<Self> {
        match b {
            b'A' | b'a' => Some(IupacRnaBase::A),
            b'C' | b'c' => Some(IupacRnaBase::C),
            b'G' | b'g' => Some(IupacRnaBase::G),
            b'U' | b'u' => Some(IupacRnaBase::U),
            b'N' | b'n' => Some(IupacRnaBase::N),
            b'R' | b'r' => Some(IupacRnaBase::R),
            b'Y' | b'y' => Some(IupacRnaBase::Y),
            b'S' | b's' => Some(IupacRnaBase::S),
            b'W' | b'w' => Some(IupacRnaBase::W),
            b'K' | b'k' => Some(IupacRnaBase::K),
            b'M' | b'm' => Some(IupacRnaBase::M),
            b'B' | b'b' => Some(IupacRnaBase::B),
            b'D' | b'd' => Some(IupacRnaBase::D),
            b'H' | b'h' => Some(IupacRnaBase::H),
            b'V' | b'v' => Some(IupacRnaBase::V),
            _ => None,
        }
    }
}

impl DnaBase {
    /// Parse a single ASCII byte into a DNA base (A,C,G,T plus IUPAC codes),
    /// case-insensitive.
    ///
    /// Returns `None` if `b` is not a valid DNA symbol.
    pub const fn from_ascii_const(b: u8) -> Option<Self> {
        match b {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'T' | b't' => Some(Self::T),
            _ => None,
        }
    }
}

impl RnaBase {
    /// Parse a single ASCII byte into a DNA base (A,C,G,T plus IUPAC codes),
    /// case-insensitive.
    ///
    /// Returns `None` if `b` is not a valid DNA symbol.
    pub const fn from_ascii_const(b: u8) -> Option<Self> {
        match b {
            b'A' | b'a' => Some(Self::A),
            b'C' | b'c' => Some(Self::C),
            b'G' | b'g' => Some(Self::G),
            b'U' | b'u' => Some(Self::U),
            _ => None,
        }
    }
}
impl Base for RnaBase {
    const ALPHABET: Alphabet = Alphabet::RNA;

    fn try_from_ascii(b: u8) -> Result<Self, BaseError> {
        if !b.is_ascii() {
            return Err(BaseError::InvalidByte {
                alphabet: Self::ALPHABET,
                invalid: b,
            });
        }

        Self::from_ascii_const(b).ok_or(BaseError::InvalidCharacter {
            alphabet: Self::ALPHABET,
            invalid: b as char,
        })
    }

    fn complement(self) -> Self {
        match self {
            Self::A => Self::U,
            Self::U => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
        }
    }

    fn to_ascii(self) -> u8 {
        match self {
            RnaBase::A => b'A',
            RnaBase::C => b'C',
            RnaBase::G => b'G',
            RnaBase::U => b'U',
        }
    }

    fn to_ascii_lower(self) -> u8 {
        match self {
            RnaBase::A => b'a',
            RnaBase::C => b'c',
            RnaBase::G => b'g',
            RnaBase::U => b'u',
        }
    }

    fn is_ambiguous(self) -> bool {
        false
    }

    fn try_chemical_class(self) -> Result<ChemClass, BaseError> {
        match self {
            RnaBase::A => Ok(ChemClass::Purine),
            RnaBase::G => Ok(ChemClass::Purine),
            RnaBase::C => Ok(ChemClass::Pyrimidine),
            RnaBase::U => Ok(ChemClass::Pyrimidine),
        }
    }

    fn to_char(self) -> char {
        self.to_ascii() as char
    }

    fn is_unambiguous(self) -> bool {
        !self.is_ambiguous()
    }
}

impl Base for IupacRnaBase {
    const ALPHABET: Alphabet = Alphabet::RNA;

    fn try_from_ascii(b: u8) -> Result<Self, BaseError> {
        if !b.is_ascii() {
            return Err(BaseError::InvalidByte {
                alphabet: Self::ALPHABET,
                invalid: b,
            });
        }

        Self::from_ascii_const(b).ok_or(BaseError::InvalidCharacter {
            alphabet: Self::ALPHABET,
            invalid: b as char,
        })
    }

    fn complement(self) -> Self {
        match self {
            Self::A => Self::U,
            Self::U => Self::A,
            Self::C => Self::G,
            Self::G => Self::C,
            Self::R => Self::Y,
            Self::Y => Self::R,
            Self::S => Self::S,
            Self::W => Self::W,
            Self::K => Self::M,
            Self::M => Self::K,
            Self::B => Self::V,
            Self::V => Self::B,
            Self::D => Self::H,
            Self::H => Self::D,
            Self::N => Self::N,
        }
    }

    fn to_ascii(self) -> u8 {
        match self {
            IupacRnaBase::A => b'A',
            IupacRnaBase::C => b'C',
            IupacRnaBase::G => b'G',
            IupacRnaBase::U => b'U',
            IupacRnaBase::N => b'N',
            IupacRnaBase::R => b'R',
            IupacRnaBase::Y => b'Y',
            IupacRnaBase::S => b'S',
            IupacRnaBase::W => b'W',
            IupacRnaBase::K => b'K',
            IupacRnaBase::M => b'M',
            IupacRnaBase::B => b'B',
            IupacRnaBase::D => b'D',
            IupacRnaBase::H => b'H',
            IupacRnaBase::V => b'V',
        }
    }

    fn to_ascii_lower(self) -> u8 {
        match self {
            IupacRnaBase::A => b'a',
            IupacRnaBase::C => b'c',
            IupacRnaBase::G => b'g',
            IupacRnaBase::U => b'u',
            IupacRnaBase::N => b'n',
            IupacRnaBase::R => b'r',
            IupacRnaBase::Y => b'y',
            IupacRnaBase::S => b's',
            IupacRnaBase::W => b'w',
            IupacRnaBase::K => b'k',
            IupacRnaBase::M => b'm',
            IupacRnaBase::B => b'b',
            IupacRnaBase::D => b'd',
            IupacRnaBase::H => b'h',
            IupacRnaBase::V => b'v',
        }
    }

    fn is_ambiguous(self) -> bool {
        !matches!(
            self,
            IupacRnaBase::A | IupacRnaBase::C | IupacRnaBase::G | IupacRnaBase::U
        )
    }

    fn try_chemical_class(self) -> Result<ChemClass, BaseError> {
        match self {
            IupacRnaBase::A => Ok(ChemClass::Purine),
            IupacRnaBase::G => Ok(ChemClass::Purine),
            IupacRnaBase::C => Ok(ChemClass::Pyrimidine),
            IupacRnaBase::U => Ok(ChemClass::Pyrimidine),
            IupacRnaBase::R => Ok(ChemClass::Purine),
            IupacRnaBase::Y => Ok(ChemClass::Pyrimidine),
            IupacRnaBase::N => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::N.to_char(),
            }),
            IupacRnaBase::S => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::S.to_char(),
            }),
            IupacRnaBase::W => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::W.to_char(),
            }),
            IupacRnaBase::K => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::K.to_char(),
            }),
            IupacRnaBase::M => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::M.to_char(),
            }),
            IupacRnaBase::B => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::B.to_char(),
            }),
            IupacRnaBase::D => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::D.to_char(),
            }),
            IupacRnaBase::H => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::H.to_char(),
            }),
            IupacRnaBase::V => Err(BaseError::AmbiguousChemicalClass {
                base: IupacRnaBase::V.to_char(),
            }),
        }
    }
}

/// Marker trait for nucleotide base types whose values are always concrete.
///
/// A `ConcreteBase` alphabet has no ambiguity symbols: every possible value of
/// the type represents exactly one biological nucleotide.
///
/// For example, `DnaBase` values are limited to `A`, `C`, `G`, and `T`, and
/// `RnaBase` values are limited to `A`, `C`, `G`, and `U`.
///
/// This trait lets ConcreteBase classes and Sequences of concrete bases provide infallible methods for
/// operations that are forced to be fallable in Base trait because of the possibility of ambiguous bases
/// such as chemical_class_concrete() or `is_palindromic() -> bool`.
///
/// This is a type-level guarantee: if `B: ConcreteBase`, then no runtime
/// ambiguity check is required.
pub trait ConcreteBase: Base + Sized {
    fn chemical_class(&self) -> ChemClass;
}

/// Marker trait for nucleotide base types whose alphabet may include
/// degenerate or ambiguous symbols.
///
/// A `DegenerateBase` alphabet can contain values that represent more than one
/// possible biological nucleotide, such as IUPAC ambiguity codes like `N`, `R`,
/// or `Y`.
///
/// Importantly, this does not mean that every value of the type is ambiguous.
/// For example, `IupacDnaBase::A` is still concrete, but `IupacDnaBase::N` is
/// ambiguous. This trait means that ambiguity is possible for this base type.
///
/// This trait lets `Seq<B>` provide fallible methods for operations that may be
/// undecidable when ambiguous bases are present, such as
/// `try_is_palindromic() -> Result<bool>`.
///
/// This is a type-level marker only. Use [`Base::is_ambiguous`] to check whether
/// a particular base value is ambiguous at runtime.
pub trait DegenerateBase: Base + Sized {
    type ConcreteEquivalent: ConcreteBase;

    fn try_to_concrete(self) -> Option<Self::ConcreteEquivalent>;
}

impl ConcreteBase for DnaBase {
    fn chemical_class(&self) -> ChemClass {
        match self {
            DnaBase::A => ChemClass::Purine,
            DnaBase::G => ChemClass::Purine,
            DnaBase::C => ChemClass::Pyrimidine,
            DnaBase::T => ChemClass::Pyrimidine,
        }
    }
}
impl ConcreteBase for RnaBase {
    fn chemical_class(&self) -> ChemClass {
        match self {
            RnaBase::A => ChemClass::Purine,
            RnaBase::G => ChemClass::Purine,
            RnaBase::U => ChemClass::Pyrimidine,
            RnaBase::C => ChemClass::Pyrimidine,
        }
    }
}
impl DegenerateBase for IupacRnaBase {
    type ConcreteEquivalent = RnaBase;

    fn try_to_concrete(self) -> Option<Self::ConcreteEquivalent> {
        match self {
            IupacRnaBase::A => Some(RnaBase::A),
            IupacRnaBase::C => Some(RnaBase::C),
            IupacRnaBase::G => Some(RnaBase::G),
            IupacRnaBase::U => Some(RnaBase::U),
            _ => None,
        }
    }
}
impl DegenerateBase for IupacDnaBase {
    type ConcreteEquivalent = DnaBase;

    fn try_to_concrete(self) -> Option<Self::ConcreteEquivalent> {
        match self {
            IupacDnaBase::A => Some(DnaBase::A),
            IupacDnaBase::C => Some(DnaBase::C),
            IupacDnaBase::G => Some(DnaBase::G),
            IupacDnaBase::T => Some(DnaBase::T),
            _ => None,
        }
    }
}

/// Const-time validator for DNA string literals.
///
/// # Panics (during const-eval)
/// Panics if `s` contains any non-DNA ASCII character.
pub const fn assert_valid_dna_literal(s: &str) {
    let bytes = s.as_bytes();
    let mut i = 0;
    while i < bytes.len() {
        if crate::base::IupacDnaBase::from_ascii_const(bytes[i]).is_none() {
            // Note: This message is intentionally simple so it works in const contexts.
            // (You can refine it later if your MSRV supports richer const panics.)
            panic!("invalid DNA base in literal");
        }
        i += 1;
    }
}

/// Names of supported alphabets.
///
/// This is mainly used for error reporting (e.g. “invalid character for DNA”).
#[derive(Debug, Clone, Copy, Eq, PartialEq)]
pub enum Alphabet {
    DNA,
    RNA,
}

impl fmt::Display for Alphabet {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Alphabet::DNA => write!(f, "DNA"),
            Alphabet::RNA => write!(f, "RNA"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dna_try_from_ascii_accepts_case_insensitive() {
        assert_eq!(IupacDnaBase::try_from_ascii(b'a').unwrap(), IupacDnaBase::A);
        assert_eq!(IupacDnaBase::try_from_ascii(b'A').unwrap(), IupacDnaBase::A);
        assert_eq!(IupacDnaBase::try_from_ascii(b't').unwrap(), IupacDnaBase::T);
    }

    #[test]
    fn rna_try_from_ascii_accepts_case_insensitive() {
        assert_eq!(IupacRnaBase::try_from_ascii(b'u').unwrap(), IupacRnaBase::U);
        assert_eq!(IupacRnaBase::try_from_ascii(b'U').unwrap(), IupacRnaBase::U);
    }

    #[test]
    fn complement_is_involutive_for_some_representative_bases() {
        // "Involutive" means comp(comp(x)) == x.
        // We only test a few representative bases to keep this minimal.
        let reps_dna = [
            IupacDnaBase::A,
            IupacDnaBase::C,
            IupacDnaBase::R,
            IupacDnaBase::B,
            IupacDnaBase::N,
        ];
        for b in reps_dna {
            assert_eq!(b.complement().complement(), b);
        }

        let reps_rna = [
            IupacRnaBase::A,
            IupacRnaBase::C,
            IupacRnaBase::R,
            IupacRnaBase::B,
            IupacRnaBase::N,
        ];
        for b in reps_rna {
            assert_eq!(b.complement().complement(), b);
        }
    }

    #[test]
    fn ascii_rendering_is_consistent() {
        // Upper then lower should match ASCII casing expectations.
        let b = IupacDnaBase::G;
        assert_eq!(b.to_ascii(), b'G');
        assert_eq!(b.to_ascii_lower(), b'g');

        let r = IupacRnaBase::U;
        assert_eq!(r.to_ascii(), b'U');
        assert_eq!(r.to_ascii_lower(), b'u');
    }

    #[test]
    fn ambiguity_flag_matches_expectations() {
        assert!(!IupacDnaBase::A.is_ambiguous());
        assert!(IupacDnaBase::N.is_ambiguous());
        assert!(IupacDnaBase::R.is_ambiguous());

        assert!(!IupacRnaBase::U.is_ambiguous());
        assert!(IupacRnaBase::N.is_ambiguous());
    }
}
