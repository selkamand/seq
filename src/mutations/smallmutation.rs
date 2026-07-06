use core::fmt;

use crate::{
    base::{Base, ChemClass, DnaBase, IupacDnaBase, IupacRnaBase, RnaBase},
    coords::{Pos, Strand},
    error::MutationError as Error,
    sequences::Seq,
};

pub(crate) type Result<T> = std::result::Result<T, Error>;

pub type IupacDnaSmallMutation = SmallMutation<IupacDnaBase>;
pub type IupacRnaSmallMutation = SmallMutation<IupacRnaBase>;
pub type DnaSmallMutation = SmallMutation<DnaBase>;
pub type RnaSmallMutation = SmallMutation<RnaBase>;

/// A small mutation (SNV/MNV/indel) over a specific nucleotide alphabet `B`.
/// - `SmallMutation<DnaBase>` is DNA
/// - `SmallMutation<RnaBase>` is RNA
///
/// ## Coordinate and allele semantics
/// - `position` is **1-based** (VCF-style) and refers to the **start** position.
/// - `reference` and `alternative` are stored **as provided** by the caller.
///   No left/right trimming, normalization, or decomposition is performed.
///
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SmallMutation<B: Base> {
    /// Chromosome / Contig of mutation
    chromosome: String,
    /// 1-based position (start)
    position: Pos,
    /// Original base sequence in reference
    reference: Seq<B>,

    /// Mutation of sequence
    alternative: Seq<B>,

    /// If mutation was on a double stranded molecule like DNA
    /// then field indicates the strand of the reference genome the reference allele was on.
    /// If mutation comes from an unstranded molecule (single stranded DNA / RNA)
    /// or strand is not known can be set to None.
    strand: Option<Strand>,

    /// - `multiallelic` indicates that this mutation originated from a site with multiple ALT
    ///   alleles (e.g. a VCF record with `ALT=A,C`). It is purely metadata for downstream
    ///   handling; it does not change coordinate or allele semantics.
    ///   It does not tell you how this was handled by mutation parsers (e.g. was alt just the
    ///   first allele, or did the parser read each allele into a separate SmallMutation object)
    multiallelic: bool,

    /// ## Filter / QC status (`pass`)
    /// - `pass` indicates whether the source record **passed upstream filtering**.
    ///   In VCF terms this typically corresponds to the `FILTER` field being `PASS`
    ///   (or `.` depending on your chosen convention).
    /// - `pass` is **metadata only**:
    ///   - it does *not* imply biological validity,
    ///   - it does *not* change classification (`SNV`/`INDEL`/etc)
    pass: bool,
}

// Implement the `fmt::Display` trait for `SmallMutation`.
impl<B: Base> fmt::Display for SmallMutation<B> {
    /// Render a compact, human-readable representation of the mutation.
    ///
    /// Format:
    /// `chrom:pos REF>ALT (delta: D; class: C; multiallelic: M; pass: P)`
    ///
    /// Intended for logging / CLI output rather than stable serialization.
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        // Use the write! macro to format the output.
        write!(
            f,
            "{}:{} {}>{} (strand: {}, delta: {}; class: {}; multiallelic:{}; pass:{})",
            self.chromosome,
            self.position,
            self.reference,
            self.alternative,
            self.strand
                .map(|opt| opt.to_string())
                .unwrap_or("None".to_string()),
            self.delta(),
            self.class(),
            self.multiallelic,
            self.pass
        )
    }
}

impl<B: Base> SmallMutation<B> {
    /// Construct a new [`SmallMutation`].
    ///
    /// This constructor does **not** attempt to validate biological correctness beyond
    /// what is already guaranteed by [`Seq<B>`] (i.e. sequences conform to the alphabet).
    ///
    /// In particular, this type:
    /// - assumes `position` is **1-based** (VCF-style coordinates)
    /// - stores `reference` and `alternative` as provided (no trimming/normalization)
    /// - treats `multiallelic` as an external flag (e.g. derived from a multi-ALT record)
    /// - treats `pass` as an external flag describing upstream filtering/QC outcome
    /// - accepts an optional `context` sequence if the caller has already computed it
    ///
    /// If you need allele normalization (left/right trimming of shared prefix/suffix),
    /// do it before constructing this type.
    ///
    /// # Parameters
    /// - `chromosome`: reference sequence / contig name (e.g. `"chr1"`)
    /// - `position`: 1-based start coordinate
    /// - `reference`: reference allele sequence
    /// - `alternative`: alternative allele sequence
    /// - `strand`: Strand
    /// - `multiallelic`: whether the originating site had multiple ALT alleles
    /// - `pass`: whether the originating record passed upstream filters/QC
    /// - `context`: optional context sequence (e.g. trinucleotide context)
    pub fn new(
        chromosome: String,
        position: Pos,
        reference: Seq<B>,
        alternative: Seq<B>,
        strand: Option<Strand>,
        multiallelic: bool,
        pass: bool,
    ) -> Self {
        Self {
            chromosome,
            position,
            reference,
            alternative,
            multiallelic,
            pass,
            strand,
        }
    }

    /// Create a minimal mutation by specifying Chromosome, Position, Reference and Alt
    /// Automatically sets other fields to a default (Strand to positive, multiallelic to false,
    /// pass to true)
    pub fn new_minimal(
        chromosome: String,
        position: Pos,
        reference: Seq<B>,
        alternative: Seq<B>,
    ) -> Self {
        Self {
            chromosome,
            position,
            reference,
            alternative,
            multiallelic: false,
            pass: true,
            strand: Some(Strand::Positive),
        }
    }

    // --- Accessors (read-only) ---

    /// Returns the chromosome / contig name (e.g. `"chr1"`).
    pub fn chromosome(&self) -> &str {
        &self.chromosome
    }

    /// Returns the 1-based start position of the mutation.
    ///
    /// This follows VCF conventions: the coordinate refers to the first base of `reference`.
    pub fn position(&self) -> Pos {
        self.position
    }

    /// Returns the reference allele sequence.
    pub fn reference(&self) -> &Seq<B> {
        &self.reference
    }

    /// Returns the alternative allele sequence.
    pub fn alternative(&self) -> &Seq<B> {
        &self.alternative
    }

    /// Returns the strand of the source molecule that contains
    /// the reference allele of this mutation (or None if from  
    /// a single stranded molecule like RNA
    pub fn strand(&self) -> Option<&Strand> {
        self.strand.as_ref()
    }

    /// Returns whether this mutation originated from a multi-allelic site.
    ///
    /// This is metadata (e.g. a VCF record with multiple ALT alleles) and does not change
    /// coordinate or allele semantics.
    pub fn is_multiallelic(&self) -> bool {
        self.multiallelic
    }

    /// Returns whether this mutation passed upstream filtering/QC.
    ///
    /// In VCF terms, this commonly corresponds to the `FILTER` field being `PASS`
    /// (and sometimes `.` depending on the caller’s convention).
    ///
    /// This is metadata only; downstream code should decide how to handle non-passing
    /// mutations explicitly.
    pub fn is_pass(&self) -> bool {
        self.pass
    }

    // --- Computed Properties (read-only) ---

    /// Return the length of the reference allele in bases.
    ///
    /// This is a convenience wrapper around [`Seq::len`].
    pub fn reflen(&self) -> usize {
        self.reference.len()
    }

    /// Return the length of the alternative allele in bases.
    ///
    /// This is a convenience wrapper around [`Seq::len`].
    pub fn altlen(&self) -> usize {
        self.alternative.len()
    }
    /// Return the signed size change implied by this mutation.
    ///
    /// Defined as:
    /// ```text
    /// delta = alt_length - ref_length
    /// ```
    ///
    /// Interpretation:
    /// - `delta == 0` → equal-length substitution (SNV/DOUBLET/MNV)
    /// - `delta > 0`  → insertion (net gain of bases)
    /// - `delta < 0`  → deletion (net loss of bases)
    ///
    /// This is purely length-based and does not depend on allele normalization.
    pub fn delta(&self) -> i64 {
        self.altlen() as i64 - self.reflen() as i64
    }
    /// Return the mutation class derived from allele lengths.
    ///
    /// This is computed on demand using [`SmallMutationType::from_lengths`]
    /// and is therefore always consistent with the stored alleles.
    pub fn class(&self) -> SmallMutationType {
        SmallMutationType::from_lengths(self.reflen(), self.altlen())
    }

    /// Compute transition/transversion classification for this mutation.
    ///
    /// This classification is only defined for **single-nucleotide substitutions**
    /// ([`SmallMutationType::SNV`]). For all other mutation types, this returns `None`.
    ///
    /// This method is conservative in the presence of ambiguity:
    /// - If either base has ambiguous chemical class (e.g. `N`, `S`, `W`), returns `None`.
    ///
    /// # Returns
    /// - `Some(TiTv::Transition)` for A↔G or C↔T substitutions (including unambiguous
    ///   IUPAC codes that resolve to a single chemical class).
    /// - `Some(TiTv::Transversion)` for purine↔pyrimidine substitutions.
    /// - `None` if not an SNV or if ambiguity prevents a confident classification.
    pub fn titv(&self) -> Option<TiTv> {
        if self.class() != SmallMutationType::SNV {
            return None;
        }

        let r = self.reference.as_slice().first()?;
        let a = self.alternative.as_slice().first()?;

        let r_chemical_class = r.try_chemical_class().ok()?;
        let a_chemical_class = a.try_chemical_class().ok()?;

        let titv = TiTv::from_chemical_class(r_chemical_class, a_chemical_class);
        Some(titv)
    }

    /// Render a compact `chrom:pos REF>ALT` representation of the mutation.
    ///
    /// This is a convenience formatter intended for **human-readable output**
    /// (e.g. logs, debugging, CLI display). It mirrors the common genomics shorthand
    /// used in VCF-style contexts, but is **not** a stable serialization format.
    ///
    /// ## Format
    ///
    /// ```text
    /// chromosome:position REF>ALT
    /// ```
    ///
    /// Where:
    /// - `chromosome` is the contig / reference name (e.g. `"chr1"`)
    /// - `position` is the **1-based** start coordinate (VCF-style)
    /// - `REF` is the reference allele sequence
    /// - `ALT` is the alternative allele sequence
    ///
    /// Alleles are rendered using their sequence display implementations
    /// (uppercase IUPAC symbols).
    ///
    /// ## Examples
    ///
    /// ```rust
    /// use seqlib::mutation::DnaSmallMutation;
    /// use seqlib::sequences::DnaSeq;
    /// use seqlib::coords::{Pos, Strand};
    ///
    /// let m = DnaSmallMutation::new(
    ///     "chr1".to_string(),
    ///     Pos::new(123).unwrap(),
    ///     DnaSeq::new("A").unwrap(),
    ///     DnaSeq::new("G").unwrap(),
    ///     Some(Strand::Positive),
    ///     false,
    ///     true,
    /// );
    ///
    /// assert_eq!(m.chrom_pos_ref_alt(), "chr1:123 A>G");
    /// ```
    pub fn chrom_pos_ref_alt(&self) -> String {
        format!(
            "{}:{} {}>{}",
            self.chromosome(),
            self.position(),
            self.reference(),
            self.alternative()
        )
    }

    /// Reverse complement mutation
    /// Reverse complements reference and alternative sequence and flips strand
    /// field if present.
    pub fn reverse_complement(&self) -> Self {
        // Flip Strand
        let newstrand = self.strand.map(|strand| strand.flip());

        Self {
            chromosome: self.chromosome.clone(),
            position: self.position,
            reference: self.reference.reverse_complement(),
            alternative: self.alternative.reverse_complement(),
            strand: newstrand,
            multiallelic: self.multiallelic,
            pass: self.pass,
        }
    }

    /// Ensure middle base of reference sequence is a pyrimidine
    /// This is accomplished by reverse complementing the mutation
    /// if there is a purine centered.
    ///
    /// # Errors
    /// If the middlebase is ambiguous, or the reference sequence is even in length (has no
    /// middlebase) will return an Error
    pub fn try_pyrimidine_center(&self) -> Result<Self> {
        let refseq = &self.reference;

        let Some(middlebase) = refseq.middlebase() else {
            return Err(Error::MissingMiddleBase);
        };

        let middle_chemclass =
            middlebase
                .try_chemical_class()
                .map_err(|_| Error::AmbiguousMiddleBase {
                    base: middlebase.to_char(),
                })?;

        // Flip Chemclass
        match middle_chemclass {
            ChemClass::Purine => Ok(self.reverse_complement()),
            ChemClass::Pyrimidine => Ok(self.clone()),
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum SmallMutationType {
    SNV,
    DOUBLET,
    MNV,
    INSERTION,
    DELETION,
}

impl fmt::Display for SmallMutationType {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            SmallMutationType::SNV => "SNV",
            SmallMutationType::DOUBLET => "DOUBLET",
            SmallMutationType::MNV => "MNV",
            SmallMutationType::INSERTION => "INSERTION",
            SmallMutationType::DELETION => "DELETION",
        };
        write!(f, "{s}")
    }
}

impl SmallMutationType {
    /// Classify a small variant based on reference and alternative allele lengths.
    ///
    /// Classification rules:
    /// - If `altlen > reflen` → [`SmallMutationType::INSERTION`]
    /// - If `altlen < reflen` → [`SmallMutationType::DELETION`]
    /// - If lengths are equal:
    ///   - `reflen == 1` → [`SmallMutationType::SNV`]
    ///   - `reflen == 2` → [`SmallMutationType::DOUBLET`]
    ///   - `reflen >= 3` → [`SmallMutationType::MNV`]
    ///
    ///
    /// # Notes
    /// - This function is purely *shape-based* and does not inspect sequence content.
    /// - A length of `0` is invalid for VCF alleles; this function currently maps
    ///   `reflen == 0 && altlen == 0` to [`SmallMutationType::MNV`] and assumes such cases are rejected
    ///   elsewhere.
    pub fn from_lengths(reflen: usize, altlen: usize) -> Self {
        match altlen.cmp(&reflen) {
            std::cmp::Ordering::Greater => Self::INSERTION,
            std::cmp::Ordering::Less => Self::DELETION,
            std::cmp::Ordering::Equal => match reflen {
                0 => Self::MNV, // OR HANDLE AS ERROR ELSEWHERE (0-LENGTH ALLELES ARE INVALID FOR VCF)
                1 => Self::SNV,
                2 => Self::DOUBLET,
                _ => Self::MNV,
            },
        }
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TiTv {
    Transition,
    Transversion,
}

impl TiTv {
    /// Classify a substitution as a transition or transversion based on chemical class.
    ///
    /// Purines -> Pyrimidine = Transversion
    /// Pyrimidine -> Purine = Transversion
    /// Purine -> Purine = Transition
    /// Pyrimidine -> Pyrimidine = Transition
    pub fn from_chemical_class(reference: ChemClass, alternative: ChemClass) -> TiTv {
        match (reference, alternative) {
            (ChemClass::Purine, ChemClass::Purine) => TiTv::Transition,
            (ChemClass::Pyrimidine, ChemClass::Pyrimidine) => TiTv::Transition,
            (ChemClass::Purine, ChemClass::Pyrimidine) => TiTv::Transversion,
            (ChemClass::Pyrimidine, ChemClass::Purine) => TiTv::Transversion,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::base::{ChemClass, IupacDnaBase, IupacRnaBase};
    use crate::sequences::Seq;

    // --- Helpers ---

    fn dna_iupac(s: &str) -> Seq<IupacDnaBase> {
        Seq::<IupacDnaBase>::new(s).unwrap()
    }

    fn _dna(s: &str) -> Seq<DnaBase> {
        Seq::<DnaBase>::new(s).unwrap()
    }

    fn rna_iupac(s: &str) -> Seq<IupacRnaBase> {
        Seq::<IupacRnaBase>::new(s).unwrap()
    }

    fn _rna(s: &str) -> Seq<RnaBase> {
        Seq::<RnaBase>::new(s).unwrap()
    }

    fn dna_mut(ref_allele: &str, alt_allele: &str) -> IupacDnaSmallMutation {
        IupacDnaSmallMutation::new(
            "chr1".to_string(),
            Pos::new(123).unwrap(),
            dna_iupac(ref_allele),
            dna_iupac(alt_allele),
            Some(Strand::Positive),
            false,
            true,
        )
    }

    fn rna_mut(ref_allele: &str, alt_allele: &str) -> IupacRnaSmallMutation {
        IupacRnaSmallMutation::new(
            "tx1".to_string(),
            Pos::new(7).unwrap(),
            rna_iupac(ref_allele),
            rna_iupac(alt_allele),
            Some(Strand::Positive),
            false,
            true,
        )
    }

    // --- Construction / context ---

    #[test]
    fn new_sets_fields_and_lengths_dna() {
        let m = dna_mut("A", "G");
        assert_eq!(m.reflen(), 1);
        assert_eq!(m.altlen(), 1);
        assert_eq!(m.delta(), 0);
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), Some(TiTv::Transition));
    }

    // --- SmallMutationType::from_lengths ---
    #[test]
    fn from_lengths_classifies_equal_length_substitutions() {
        assert_eq!(
            SmallMutationType::from_lengths(1, 1),
            SmallMutationType::SNV
        );
        assert_eq!(
            SmallMutationType::from_lengths(2, 2),
            SmallMutationType::DOUBLET
        );
        assert_eq!(
            SmallMutationType::from_lengths(3, 3),
            SmallMutationType::MNV
        );
        assert_eq!(
            SmallMutationType::from_lengths(10, 10),
            SmallMutationType::MNV
        );
    }

    #[test]
    fn from_lengths_classifies_indels() {
        assert_eq!(
            SmallMutationType::from_lengths(1, 2),
            SmallMutationType::INSERTION
        );
        assert_eq!(
            SmallMutationType::from_lengths(2, 1),
            SmallMutationType::DELETION
        );
        assert_eq!(
            SmallMutationType::from_lengths(5, 9),
            SmallMutationType::INSERTION
        );
        assert_eq!(
            SmallMutationType::from_lengths(9, 5),
            SmallMutationType::DELETION
        );
    }

    // --- delta / class integration tests ---

    #[test]
    fn class_and_delta_match_expected_for_snv_doublet_mnv() {
        let snv = dna_mut("A", "C");
        assert_eq!(snv.class(), SmallMutationType::SNV);
        assert_eq!(snv.delta(), 0);

        let dbl = dna_mut("AC", "GT");
        assert_eq!(dbl.class(), SmallMutationType::DOUBLET);
        assert_eq!(dbl.delta(), 0);

        let mnv = dna_mut("ACG", "TTA");
        assert_eq!(mnv.class(), SmallMutationType::MNV);
        assert_eq!(mnv.delta(), 0);
    }

    #[test]
    fn class_and_delta_match_expected_for_insertion_and_deletion() {
        let ins = dna_mut("A", "AT");
        assert_eq!(ins.class(), SmallMutationType::INSERTION);
        assert_eq!(ins.delta(), 1);

        let del = dna_mut("AT", "A");
        assert_eq!(del.class(), SmallMutationType::DELETION);
        assert_eq!(del.delta(), -1);
    }

    // --- TiTv::from_chemical_class ---

    #[test]
    fn titv_from_chemical_class_transition_and_transversion() {
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Purine, ChemClass::Purine),
            TiTv::Transition
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Pyrimidine, ChemClass::Pyrimidine),
            TiTv::Transition
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Purine, ChemClass::Pyrimidine),
            TiTv::Transversion
        );
        assert_eq!(
            TiTv::from_chemical_class(ChemClass::Pyrimidine, ChemClass::Purine),
            TiTv::Transversion
        );
    }

    // --- SmallMutation::titv ---

    #[test]
    fn titv_only_defined_for_snvs() {
        let mnv = dna_mut("AC", "GT");
        assert_eq!(mnv.class(), SmallMutationType::DOUBLET);
        assert_eq!(mnv.titv(), None);

        let ins = dna_mut("A", "AT");
        assert_eq!(ins.class(), SmallMutationType::INSERTION);
        assert_eq!(ins.titv(), None);

        let del = dna_mut("AT", "A");
        assert_eq!(del.class(), SmallMutationType::DELETION);
        assert_eq!(del.titv(), None);
    }

    #[test]
    fn titv_transition_examples_dna() {
        // A <-> G is a transition (purine <-> purine)
        let m = dna_mut("A", "G");
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), Some(TiTv::Transition));

        // C <-> T is a transition (pyrimidine <-> pyrimidine)
        let m2 = dna_mut("C", "T");
        assert_eq!(m2.titv(), Some(TiTv::Transition));
    }

    #[test]
    fn titv_transversion_examples_dna() {
        // A <-> C is a transversion (purine <-> pyrimidine)
        let m = dna_mut("A", "C");
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), Some(TiTv::Transversion));

        // G <-> T is a transversion
        let m2 = dna_mut("G", "T");
        assert_eq!(m2.titv(), Some(TiTv::Transversion));
    }

    #[test]
    fn titv_returns_none_when_ambiguous_base_present() {
        // N is ChemClass::Ambiguous in your design
        let m = dna_mut("N", "A");
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), None);

        let m2 = dna_mut("A", "N");
        assert_eq!(m2.class(), SmallMutationType::SNV);
        assert_eq!(m2.titv(), None);
    }

    // --- Display ---

    #[test]
    fn display_includes_core_fields() {
        let m = dna_mut("A", "G");
        let s = m.to_string();

        // Keep this intentionally loose so formatting tweaks don't require rewrites.
        assert!(s.contains("chr1:123"));
        assert!(s.contains("A>G"));
        assert!(s.contains("delta: 0"));
        assert!(s.contains("class: SNV"));
        assert!(s.contains("multiallelic:false"));
    }

    // --- RNA smoke tests (generic over Base works) ---

    #[test]
    fn rna_small_mutation_works_and_classifies() {
        let m = rna_mut("A", "G");
        assert_eq!(m.reflen(), 1);
        assert_eq!(m.altlen(), 1);
        assert_eq!(m.delta(), 0);
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), Some(TiTv::Transition)); // purine<->purine is transition

        let ins = rna_mut("A", "AU");
        assert_eq!(ins.class(), SmallMutationType::INSERTION);
        assert_eq!(ins.delta(), 1);
        assert_eq!(ins.titv(), None);
    }

    #[test]
    fn rna_ambiguous_titv_none() {
        // 'N' exists for RNA alphabet in your base impl
        let m = rna_mut("N", "A");
        assert_eq!(m.class(), SmallMutationType::SNV);
        assert_eq!(m.titv(), None);
    }
}
