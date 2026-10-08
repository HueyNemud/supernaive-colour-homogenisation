//! Project settings: `homog.toml`, every key optional and overridden by the command line.
//!
//! ```toml
//! target = [255, 255, 255]
//! lighting = "normal"          # even | normal | uneven | none
//! icc = true
//! [hints]
//! paper = "hints/paper"        # imagettes of bare paper
//! keep = "hints/keep"          # imagettes of colours that must not become white
//! [seams]
//! gcp_dir = "GCP"              # QGIS .points files giving the layout of the sheets
//! [output]
//! dir = "output"
//! format = "auto"              # auto | tif | png | jpg
//! jpeg_quality = 95
//! compress = "DEFLATE"
//! ```
//! Relative paths are relative to the configuration file.

use std::path::{Path, PathBuf};

use anyhow::{Context, Result};
use serde::Deserialize;

use crate::background::Lighting;
use crate::io::Format;

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    pub target: Option<[u8; 3]>,
    pub lighting: Option<Lighting>,
    pub icc: Option<bool>,
    #[serde(default)]
    pub hints: HintsConfig,
    #[serde(default)]
    pub seams: SeamsConfig,
    #[serde(default)]
    pub output: OutputConfig,
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct HintsConfig {
    pub paper: Option<PathBuf>,
    pub keep: Option<PathBuf>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SeamsConfig {
    pub gcp_dir: Option<PathBuf>,
}

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct OutputConfig {
    pub dir: Option<PathBuf>,
    pub format: Option<Format>,
    pub jpeg_quality: Option<u8>,
    pub compress: Option<String>,
}

impl Config {
    pub fn load(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).with_context(|| format!("cannot read {}", path.display()))?;
        Self::parse(&text, path.parent().unwrap_or(Path::new(""))).with_context(|| format!("invalid configuration {}", path.display()))
    }

    /// Parse a configuration, relative paths being resolved against `base`.
    pub fn parse(text: &str, base: &Path) -> Result<Self> {
        let mut c: Config = toml::from_str(text)?;
        for p in [&mut c.hints.paper, &mut c.hints.keep, &mut c.seams.gcp_dir, &mut c.output.dir].into_iter().flatten() {
            if p.is_relative() {
                *p = base.join(&*p);
            }
        }
        Ok(c)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_and_resolves_paths() {
        let c = Config::parse(
            "lighting = \"uneven\"\n[hints]\nkeep = \"hints/keep\"\n[seams]\ngcp_dir = \"/abs/GCP\"\n[output]\nformat = \"png\"\n",
            Path::new("/project"),
        )
        .unwrap();
        assert_eq!(c.lighting, Some(Lighting::Uneven));
        assert_eq!(c.hints.keep.as_deref(), Some(Path::new("/project/hints/keep")));
        assert_eq!(c.hints.paper, None);
        assert_eq!(c.seams.gcp_dir.as_deref(), Some(Path::new("/abs/GCP")));
        assert_eq!(c.output.format, Some(Format::Png));
    }

    #[test]
    fn rejects_unknown_keys() {
        assert!(Config::parse("degree = 3\n", Path::new("")).is_err());
        assert!(Config::parse("[hints]\npapers = \"x\"\n", Path::new("")).is_err());
    }
}
