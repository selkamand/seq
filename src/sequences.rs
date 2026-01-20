use core::fmt;

/// Seq class for sequence representations
#[derive(Debug)]
pub struct Seq {
    seq: String,
    alphabet: Alphabet,
}

impl fmt::Display for Seq {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{}", self.seq)
    }
}

impl Seq {
    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn is_empty(&self) -> bool {
        self.seq.is_empty()
    }

    pub fn new(seq: String, alphabet: Alphabet) -> Result<Self, &'static str> {
        Ok(Self { seq, alphabet })
    }
}

#[derive(Debug)]
enum Alphabet {
    DNA,
    RNA,
}

impl Alphabet {
    fn is_valid_char(&self, c: char) -> bool {
        match self {
            Alphabet::DNA => matches!(c, 'A' | 'C' | 'G' | 'T' | 'N'),
            Alphabet::RNA => matches!(c, 'A' | 'C' | 'G' | 'U' | 'N'),
        }
    }
}
