use crate::{
    base::Base,
    coords::{Interval, Pos},
    error::MutationError,
    mutations::SmallMutation,
    sequences::{Seq, SourcedSeq},
};

pub(crate) type Result<T> = std::result::Result<T, MutationError>;

/// A small mutation annotated with the reference sequence from which it was created
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct MutationWithContext<B: Base> {
    mutation: SmallMutation<B>,
    context: SourcedSeq<B>,
    anchor: Pos,
}

impl<B: Base> std::fmt::Display for MutationWithContext<B> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        writeln!(f, "-- Mutation --")?;
        writeln!(f, "{}", self.mutation)?;
        writeln!(f)?;
        writeln!(f, "-- Context --")?;
        writeln!(f, "{}", self.context)
    }
}

impl<B: Base> MutationWithContext<B> {
    /// Construct a new [`MutationWithContext`] struct
    ///
    /// ## Errors
    /// May return a [`MutationError::MismatchedStrandOption`] error if a strand is specified
    /// for mutation or context but not both.
    ///
    /// Note if strand is supplied for both but the value itself differs, the context will be
    /// reverse complemented so they match
    pub fn new(mutation: SmallMutation<B>, context: SourcedSeq<B>, anchor: Pos) -> Result<Self> {
        // Strand comparison
        let mutstrand = mutation.strand();
        let contextstrand = context.strand();

        //TODO: Because this function takes ownership of context, we could use
        //reverse_complement_in_place to make more memory efficient but because
        //contexts will usually be quite short and very few in memory at any one time we will
        //keep the copy-on-modify approach here for now
        let normalised_context = match (mutstrand, contextstrand) {
            (None, None) => context,

            // If strand is given for mutation but not for context sequence or vice-versa throw an error
            (None, Some(_)) => return Err(MutationError::MismatchedStrandOption),
            (Some(_), None) => return Err(MutationError::MismatchedStrandOption),

            // Reverse complement context to ensure strand of context and mutation always match
            (Some(strand1), Some(strand2)) => match strand1 == strand2 {
                true => context,
                false => context.reverse_complement(),
            },
        };

        Ok(Self {
            mutation,
            context: normalised_context,
            anchor,
        })
    }

    // < Getters >

    /// Sequence-local 1-based position within `context sequence` describing where the mutation
    /// starts.
    ///
    /// For example, `anchor of Pos(1) means mutation of context() `seq[0]`.
    ///
    pub fn anchor(&self) -> &Pos {
        &self.anchor
    }

    /// Return the anchor index as a 0-based index into `seq`.
    ///
    /// Returns `None` if the anchor lies outside the stored sequence.
    pub fn anchor_index0(&self) -> Option<usize> {
        let idx0 = self.anchor().get().checked_sub(1)?;
        (idx0 < self.context().seq().len()).then_some(idx0)
    }

    /// Get the mutation
    pub fn mutation(&self) -> &SmallMutation<B> {
        &self.mutation
    }

    /// Get the context of the mutation
    pub fn context(&self) -> &SourcedSeq<B> {
        &self.context
    }

    // pub fn ref_trinuc(&self) -> Option<&[B]> {
    //     self.context.as_ref()?.kmer_centered_on_anchor(3)
    // }
    //
    // pub fn ref_pentanuc(&self) -> Option<&[B]> {
    //     if self.mutation.class() != SmallMutationType::SNV {
    //         return None;
    //     }
    //     self.context.as_ref()?.kmer_centered_on_anchor(5)
    // }

    /// Get the interval of context affected by the mutation
    pub fn mutated_interval(&self) -> Result<Interval> {
        let pos1 = self.anchor().to_owned();
        let pos2 = pos1.try_add(self.mutation().reflen().saturating_sub(1))?;

        let interval = Interval::new(pos1, pos2)?;

        Ok(interval)
    }

    /// Apply a mutation to a sequence context
    pub fn apply_mutation(&self) -> Result<Seq<B>> {
        let mut seq = self.context().seq().clone();

        let interval = self.mutated_interval()?;

        seq.mutate(interval, self.mutation().alternative())
            .map_err(|source| MutationError::FailedToApplyMutationToContext {
                mutation: self.mutation().chrom_pos_ref_alt(),
                source,
            })
    }

    /// Borrow a reference k-mer of length `k` centered on the anchor position.
    ///
    /// The anchor base lies exactly at the center of the returned k-mer.
    /// `k` must be **odd and non-zero**.
    ///
    /// ## Returns
    /// - `Some(&[B])` if a centered k-mer of length `k` can be represented
    ///   by the stored context window
    /// - `None` otherwise
    ///
    /// ## `None` is returned when:
    /// - `k` is zero or even
    /// - the anchor lies outside the stored sequence
    /// - there is insufficient flanking sequence on either side of the anchor
    ///
    /// ## Notes
    /// - No mutation-type checks are performed (e.g. SNV vs indel)
    /// - No strand normalization or reverse-complementing is applied
    /// - The returned slice is borrowed; no allocation or copying occurs
    pub fn kmer_centered_on_anchor(&self, k: usize) -> Option<&[B]> {
        if k == 0 || k.is_multiple_of(2) {
            return None;
        }
        let center = self.anchor_index0()?;
        let half = k / 2;
        let start = center.checked_sub(half)?;
        let end = center + half + 1;
        if end > self.context().seq().len() {
            return None;
        }
        // subseq_slice returns Result<&[B]>; convert to Option
        self.context().seq().slice(start, end).ok()
    }

    /// Visualise the mutation in context
    ///
    /// creates a string with two lines. One with reference sequence (brakets indicating mutated
    /// interval/position) and one with the alternative sequence
    ///
    /// If there were errors applying the mutation to the sequence context, alternative will just
    /// display the error message.
    pub fn to_difference_string(&self) -> String {
        // Grab context
        let ctx = self.context();

        // Create formatted reference string
        let refstring = if let Ok(interval) = self.mutated_interval() {
            ctx.seq().format_with_highlight_interval(Some(&interval))
        } else {
            ctx.seq()
                .format_with_highlight_pos(Some(self.mutation().position()))
        };

        // Create formatted altstring
        let mutated_sequence = self.apply_mutation();

        let altstring = match mutated_sequence {
            Ok(alt) => alt.to_string(),
            Err(e) => format!("Failed to apply mutation: {e}"),
        };

        format!("{refstring}\n{altstring}")
    }
}
