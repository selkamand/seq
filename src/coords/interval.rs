//! Interval Type (simple range defined by start and end position)

use crate::coords::Pos;
use crate::error::CoordError as Error;
pub(crate) type Result<T> = std::result::Result<T, Error>;

/// A genomic interval (Start & End)
/// Both are 1-based and both-end inclusive
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Interval {
    start: Pos,
    end: Pos,
}

impl std::fmt::Display for Interval {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}-{}", self.start, self.end)
    }
}

impl Interval {
    /// Creates a 1-based inclusive interval from `start` to `end`.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{Interval, Pos};
    ///
    /// let interval = Interval::new(Pos::new(2)?, Pos::new(5)?)?;
    ///
    /// assert_eq!(*interval.start(), Pos::new(2)?);
    /// assert_eq!(*interval.end(), Pos::new(5)?);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    ///
    /// # Errors
    ///
    /// Returns [`Error::RangeEndTooSmall`] if `end` is less than `start`.
    pub fn new(start: Pos, end: Pos) -> Result<Self> {
        if end < start {
            return Err(Error::RangeEndTooSmall { start, end });
        }
        Ok(Self { start, end })
    }

    /// Creates an interval around `pos` with `left` bases before and `right` bases after it.
    ///
    /// The bounds saturate at [`Pos::MIN`] and [`Pos::MAX`] rather than failing on
    /// underflow or overflow.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{Interval, Pos};
    ///
    /// let interval = Interval::around_position(Pos::new(10)?, 2, 3);
    ///
    /// assert_eq!(*interval.start(), Pos::new(8)?);
    /// assert_eq!(*interval.end(), Pos::new(13)?);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn around_position(pos: Pos, left: usize, right: usize) -> Self {
        Self {
            start: pos.saturating_sub(left),
            end: pos.saturating_add(right),
        }
    }

    /// Returns the inclusive start position.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{Interval, Pos};
    ///
    /// let interval = Interval::new(Pos::new(3)?, Pos::new(7)?)?;
    ///
    /// assert_eq!(*interval.start(), Pos::new(3)?);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn start(&self) -> &Pos {
        &self.start
    }

    /// Returns the inclusive end position.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{Interval, Pos};
    ///
    /// let interval = Interval::new(Pos::new(3)?, Pos::new(7)?)?;
    ///
    /// assert_eq!(*interval.end(), Pos::new(7)?);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn end(&self) -> &Pos {
        &self.end
    }

    /// Check if region is empty. Always returns false as regions are never empty, by definition they contain at least 1 base)
    pub fn is_empty(&self) -> bool {
        false
    }

    /// Returns the number of positions spanned by the interval.
    ///
    /// Because intervals are inclusive, `1-1` has length 1.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{Interval, Pos};
    ///
    /// let interval = Interval::new(Pos::new(2)?, Pos::new(5)?)?;
    ///
    /// assert_eq!(interval.len(), 4);
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
    pub fn len(&self) -> usize {
        self.end
            .get()
            .saturating_sub(self.start.get())
            .saturating_add(1)
    }

    /// Returns the interval as 0-based half-open indices.
    ///
    /// The returned `(start, end)` pair is suitable for Rust slicing, where
    /// `start` is inclusive and `end` is exclusive.
    ///
    /// # Examples
    ///
    /// ```
    /// use seqlib::coords::{Interval, Pos};
    ///
    /// let interval = Interval::new(Pos::new(2)?, Pos::new(5)?)?;
    ///
    /// assert_eq!(interval.as_0based_indices(), (1, 5));
    /// # Ok::<(), seqlib::error::CoordError>(())
    /// ```
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
