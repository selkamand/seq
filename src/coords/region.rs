//!  Region type

use crate::coords::Interval;

/// A region (name of originating sequence + 1-based inclusive interval)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    pub name: String,
    pub interval: Interval,
}

impl std::fmt::Display for Region {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}:{}-{}",
            self.name, self.interval.start, self.interval.end
        )
    }
}

impl Region {
    /// Create a new [`Region`] from a sequence **name** and [`Interval`]
    pub fn new<N, I>(name: N, interval: I) -> Self
    where
        N: Into<String>,
        I: Into<Interval>,
    {
        Self {
            name: name.into(),
            interval: interval.into(),
        }
    }
}
