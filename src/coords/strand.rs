//! Strand definition

/// The strand. Exactly what this enum means
/// can depend on the context.
#[derive(Debug, Clone, PartialEq, Eq, Copy)]
pub enum Strand {
    /// The positive strand (`+`).
    ///
    /// In some contexts called the _sense_ or _forward_ strand.
    Positive,
    /// The negative strand (`-`).
    ///
    /// In some contexts described as the _antisense_ or _reverse_ strand.
    Negative,
}

impl std::fmt::Display for Strand {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Strand::Positive => write!(f, "+"),
            Strand::Negative => write!(f, "-"),
        }
    }
}

impl Strand {
    /// Flip strand from Positive to Negative or Negative to Positive
    pub fn flip(&self) -> Self {
        match self {
            Strand::Positive => Strand::Negative,
            Strand::Negative => Strand::Positive,
        }
    }
}
