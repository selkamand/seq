//! # seqlib
//!
//! `seqlib` is a small, type-safe Rust library for working with DNA and RNA sequences.
//!
//! It provides a robust core representation for biological sequences with:
//! - explicit DNA / RNA alphabets
//! - full IUPAC ambiguity support
//! - conservative, correctness-first semantics
//!
//! ## Design philosophy
//!
//! - correctness over convenience
//! - explicit handling of ambiguous bases
//! - no I/O, no parsing, no side effects
//! - intended as a foundation for larger bioinformatics pipelines
//!
//! ## Module overview
//!
//! - [`sequences`] — validated DNA/RNA sequences and core transformations
//! - [`base`] — nucleotide alphabets
//! - [`context`] — reference sequence windows around mutations
//! - [`mutation`] — small variants and mutation metadata
//! - [`coords`] — safe 1-based biological coordinates
//!
//! ## What this crate does *not* do
//!
//! - FASTA/FASTQ parsing
//! - file I/O
//! - alignment or variant calling
//!
//! Those concerns are intentionally left to downstream crates.

pub mod base;
pub mod context;
pub mod coords;
pub mod error;
pub mod mutations;
pub mod sequences;
