//!  Region type

use crate::coords::Interval;

/// A region (name of originating sequence + 1-based inclusive interval)
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Region {
    name: String,
    interval: Interval,
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
