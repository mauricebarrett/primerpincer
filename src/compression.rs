use clap::ValueEnum;
use niffler::Level;
use niffler::send::compression::Format;

/// Output compression formats supported by PrimerPincer.
/// These cover the common bioinformatics compression formats supported by niffler.
#[derive(Copy, Clone, Debug, Eq, PartialEq, ValueEnum)]
pub enum OutputCompression {
    /// No compression; write plain text FASTQ.
    None,
    /// Standard gzip compression.
    Gzip,
    /// bzip2 compression.
    #[value(alias = "bz2", alias = "bzip")]
    Bzip2,
    /// LZMA/XZ compression.
    #[value(alias = "lzma")]
    Xz,
    /// Zstandard compression.
    #[value(alias = "zst")]
    Zstd,
}

impl OutputCompression {
    /// Map the enum to the corresponding niffler compression format.
    pub fn to_format(self) -> Format {
        match self {
            OutputCompression::None => Format::No,
            OutputCompression::Gzip => Format::Gzip,
            OutputCompression::Bzip2 => Format::Bzip,
            OutputCompression::Xz => Format::Lzma,
            OutputCompression::Zstd => Format::Zstd,
        }
    }

    /// Choose a reasonable default compression level for the format.
    pub fn default_level(self) -> Level {
        match self {
            OutputCompression::None => Level::One,
            // Level::One is a good default for gzip
            OutputCompression::Gzip => Level::One,
            // Level::Nine for stronger bzip2 compression
            OutputCompression::Bzip2 => Level::Nine,
            // Level::Six is a balanced default for xz/lzma and zstd
            OutputCompression::Xz => Level::Six,
            OutputCompression::Zstd => Level::Six,
        }
    }
}
