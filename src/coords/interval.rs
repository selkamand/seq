//! Interval Type (simple range defined by start and end position)

use crate::coords::Pos;
use crate::error::CoordError as Error;
pub(crate) type Result<T> = std::result::Result<T, Error>;

/// A genomic interval (Start & End)
/// Both are 1-based and both-end inclusive
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    pub start: Pos,
    pub end: Pos,
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", self.start, self.end)
    }
}

impl Interval {
    /// Create a new region (1-based inclusive)
    pub fn new(start: Pos, end: Pos) -> Result<Self> {
        if end < start {
            return Err(Error::RangeEndTooSmall { start, end });
        }
        Ok(Self { start, end })
    }

    /// Get start position
    pub fn start(&self) -> Pos {
        self.start
    }

    /// Get end position
    pub fn end(&self) -> Pos {
        self.end
    }

    /// Check if region is empty. Always returns false as regions are never empty, by definition they contain at least 1 base)
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Length of region (number of bases they span)
    pub fn len(&self) -> usize {
        self.end
            .get()
            .saturating_sub(self.start.get())
            .saturating_add(1)
    }

    /// Get region as 0-based indices for easier extractions of sequence slices
    pub fn as_0based_indices(&self) -> (usize, usize) {
        (self.start.as_0based_index(), self.end.as_0based_index() + 1)
    }
}

impl Default for Interval {
    fn default() -> Self {
        Self {
            start: Pos::MIN,
            end: Pos::MIN,
        }
    }
}
