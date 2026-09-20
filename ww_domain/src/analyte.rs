//! Multi-target analyte taxonomy for wastewater-based epidemiology (WBE).
//!
//! The physical transport equation governing analyte fate in the sewer network
//! is identical across all surveillance targets — only the **first-order decay
//! rate** `k` (d⁻¹) and the **source generation term** `f` differ between
//! analyte classes.  This module encodes that taxonomy so every downstream
//! component (EWMA α selection, detection threshold, severity policy, SHPINN
//! regression pipelines) can adapt without hard-coding analyte names.
//!
//! ```text
//! SHPINN transport (1-D advection–diffusion–decay):
//!
//!   ∂C/∂t + v ∂C/∂x = D ∂²C/∂x² − k·C + f(x,t)
//!
//!   C : analyte concentration, LINEAR scale [copies/L or µg/L]
//!       (the PDE is linear in C; for ℓ = log₁₀ C it becomes
//!        ℓₜ + vℓₓ = Dℓₓₓ + D ln10 ℓₓ² − k/ln10 + f/(C ln10))
//!   v : mean flow velocity
//!   D : effective axial dispersion
//!   k : first-order decay rate constant (THIS file encodes k per analyte)
//!   f : source generation term       (encoded as baseline_log + outbreak pulse)
//! ```
//!
//! ## Multi-target catalog
//!
//! | Category               | Example analytes           | k range (d⁻¹) | Primary application |
//! |------------------------|----------------------------|---------------|---------------------|
//! | InfectiousPathogen     | SARS-CoV-2, Influenza, …   | 0.30–0.55     | Outbreak early-warning |
//! | PharmaceuticalAMR      | Amoxicillin, Ciprofloxacin | 0.05–0.20     | AMR hotspot mapping |
//! | IllicitSubstance       | Fentanyl, Cocaine, MDMA    | 0.20–0.28     | Inverse source localisation |
//! | IndustrialToxicant     | PFAS, Pb, TCE              | 0.002–0.04    | Illegal dumping tracing |

use serde::{Deserialize, Serialize};

// ── Category enum ─────────────────────────────────────────────────────────────

/// Surveillance target category — governs kinetic parameters, EWMA α
/// selection, detection thresholds, and severity policy.
#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub enum AnalyteCategory {
    /// Viruses, bacteria, and parasites with temperature-dependent viral
    /// or microbial decay kinetics (`k_pathogen`, typically 0.30–0.55 d⁻¹).
    /// Primary application: outbreak early-warning and community prevalence tracking.
    InfectiousPathogen,

    /// Prescription pharmaceuticals and antimicrobial resistance (AMR) markers
    /// governed by chemical half-life and hydrolysis kinetics
    /// (`k_chem`, typically 0.05–0.20 d⁻¹).
    /// Primary application: AMR hotspot mapping and over-prescription surveillance.
    PharmaceuticalAMR,

    /// Illicit drugs and their urinary/biliary metabolites, characterised by
    /// high-throughput LC-MS/MS quantitation of stable metabolic products
    /// (`k_metabolite`, typically 0.20–0.28 d⁻¹).
    /// Primary application: inverse source localisation for public-health / law
    /// enforcement interventions.
    IllicitSubstance,

    /// Heavy metals, PFAS, and industrial solvents — near-zero or very slow
    /// decay (k ≈ 0–0.04 d⁻¹) due to environmental persistence or sorption
    /// to biofilm / suspended solids.
    /// Primary application: illegal industrial dumping tracing.
    IndustrialToxicant,
}

impl AnalyteCategory {
    /// Lowercase snake-case identifier — stored in the ontology as a string
    /// property so it survives serialisation round-trips without depending on
    /// the Rust enum's discriminant ordering.
    pub fn as_str(&self) -> &'static str {
        match self {
            AnalyteCategory::InfectiousPathogen => "infectious_pathogen",
            AnalyteCategory::PharmaceuticalAMR  => "pharmaceutical_amr",
            AnalyteCategory::IllicitSubstance   => "illicit_substance",
            AnalyteCategory::IndustrialToxicant => "industrial_toxicant",
        }
    }

    pub fn from_str(s: &str) -> Self {
        match s {
            "pharmaceutical_amr"  => AnalyteCategory::PharmaceuticalAMR,
            "illicit_substance"   => AnalyteCategory::IllicitSubstance,
            "industrial_toxicant" => AnalyteCategory::IndustrialToxicant,
            _                     => AnalyteCategory::InfectiousPathogen,
        }
    }

    /// Human-readable label for reports and terminal output.
    pub fn display_label(&self) -> &'static str {
        match self {
            AnalyteCategory::InfectiousPathogen => "Infectious Pathogen",
            AnalyteCategory::PharmaceuticalAMR  => "Pharmaceutical / AMR",
            AnalyteCategory::IllicitSubstance   => "Illicit Substance",
            AnalyteCategory::IndustrialToxicant => "Industrial Toxicant",
        }
    }

    /// Terminal colour prefix for category-differentiated output.
    pub fn tag(&self) -> &'static str {
        match self {
            AnalyteCategory::InfectiousPathogen => "[PATH]",
            AnalyteCategory::PharmaceuticalAMR  => "[PHRM]",
            AnalyteCategory::IllicitSubstance   => "[ILCT]",
            AnalyteCategory::IndustrialToxicant => "[INDX]",
        }
    }
}

// ── AnalyteProfile ────────────────────────────────────────────────────────────

/// Complete kinetic + detection profile for one surveillance analyte.
///
/// The `decay_rate_k` field feeds two consumers:
///
/// 1. **EWMA α** — derived as `α = clamp(1 − e^{−k}, 0.10, 0.40)`.  Fast-
///    decaying analytes (high k) get high α so the baseline reacts quickly;
///    persistent analytes (k ≈ 0) are clamped to 0.10 for stability.
///
/// 2. **SHPINN regression pipelines** — exposed through the ontology
///    `decay_rate_k` property on every `PathogenSignal` so downstream
///    PDE solvers can read the correct k without hard-coding.
#[derive(Debug, Clone)]
pub struct AnalyteProfile {
    /// Human-readable analyte name (e.g., `"SARS-CoV-2"`, `"Fentanyl"`).
    pub name: &'static str,

    /// Detection target identifier: gene (qPCR), compound (LC-MS/MS),
    /// isotope (ICP-MS), etc.
    pub target_marker: &'static str,

    /// Analytical method string (e.g., `"qPCR"`, `"LC-MS/MS"`, `"ICP-MS"`).
    pub method: &'static str,

    /// Surveillance category — governs kinetics and detection policy.
    pub category: AnalyteCategory,

    /// Endemic background concentration in log₁₀ (copies/L or µg/L),
    /// site-independent.  Outbreak pulses are added on top of this baseline.
    pub baseline_log: f64,

    /// Day-to-day noise standard deviation in log₁₀ space.
    pub noise_std: f64,

    /// First-order decay rate constant in wastewater (d⁻¹).
    /// Near-zero for persistent contaminants (PFAS, heavy metals).
    /// This is the `k` parameter in the SHPINN transport PDE.
    pub decay_rate_k: f64,

    /// Minimum z-score to raise any alert for this analyte.
    ///
    /// Pathogens use 2.5 (early-warning priority).
    /// Stable categories use higher thresholds to reduce false positives.
    pub z_threshold: f64,
}

impl AnalyteProfile {
    /// Derive the EWMA smoothing factor α from the analyte's decay rate.
    ///
    /// Derivation: in one observation period (1 day) a fraction `1 − e^{−k}`
    /// of the signal decays, so the baseline must have updated by at least that
    /// fraction to remain representative of the current true level.
    ///
    /// ```text
    /// α = clamp(1 − e^{−k}, 0.10, 0.40)
    /// ```
    ///
    /// | Category            | k (d⁻¹) | raw α | clamped α |
    /// |---------------------|----------|-------|-----------|
    /// | Infectious pathogen | 0.50     | 0.39  | 0.39      |
    /// | Illicit substance   | 0.25     | 0.22  | 0.22      |
    /// | Pharmaceutical      | 0.10     | 0.095 | 0.10      |
    /// | PFAS (industrial)   | 0.005    | 0.005 | 0.10      |
    pub fn ewma_alpha(&self) -> f64 {
        let raw = 1.0_f64 - (-self.decay_rate_k).exp();
        raw.clamp(0.10, 0.40)
    }
}

// ── Multi-target analyte catalog ──────────────────────────────────────────────
//
// To add a new analyte: append an `AnalyteProfile` entry here.  All
// downstream systems (scenario generator, detector, schema) pick it up
// automatically.

/// Canonical analyte catalog for the wastewater biosurveillance system.
///
/// The 18 analytes span all four `AnalyteCategory` classes and represent
/// the four primary application domains described in the SHPINN specification.
/// Add or remove entries here to customise the surveillance panel.
pub const ANALYTE_CATALOG: &[AnalyteProfile] = &[

    // ═══════════════════════════════════════════════════════════════════════
    // INFECTIOUS PATHOGENS — temperature-dependent viral / bacterial decay
    // k_pathogen: 0.30–0.55 d⁻¹   α: 0.26–0.39
    // ═══════════════════════════════════════════════════════════════════════

    AnalyteProfile {
        name: "SARS-CoV-2", target_marker: "N1", method: "qPCR",
        category:     AnalyteCategory::InfectiousPathogen,
        baseline_log: 3.8,  noise_std: 0.22,
        decay_rate_k: 0.50, z_threshold: 2.5,
    },
    AnalyteProfile {
        name: "Mpox", target_marker: "E6L", method: "ddPCR",
        category:     AnalyteCategory::InfectiousPathogen,
        baseline_log: 1.0,  noise_std: 0.18,
        decay_rate_k: 0.45, z_threshold: 2.5,
    },
    AnalyteProfile {
        name: "Vibrio_cholerae", target_marker: "ctxA", method: "qPCR",
        category:     AnalyteCategory::InfectiousPathogen,
        baseline_log: 2.2,  noise_std: 0.20,
        decay_rate_k: 0.40, z_threshold: 2.5,
    },
    AnalyteProfile {
        name: "Influenza_A", target_marker: "M-gene", method: "RT-qPCR",
        category:     AnalyteCategory::InfectiousPathogen,
        baseline_log: 2.5,  noise_std: 0.25,
        decay_rate_k: 0.55, z_threshold: 2.5,
    },
    AnalyteProfile {
        name: "Poliovirus", target_marker: "5-UTR", method: "RT-qPCR",
        category:     AnalyteCategory::InfectiousPathogen,
        baseline_log: 0.5,  noise_std: 0.15,
        decay_rate_k: 0.35, z_threshold: 2.5,
    },
    AnalyteProfile {
        name: "Norovirus_GII", target_marker: "ORF1-ORF2", method: "RT-ddPCR",
        category:     AnalyteCategory::InfectiousPathogen,
        baseline_log: 3.2,  noise_std: 0.28,
        decay_rate_k: 0.30, z_threshold: 2.5,
    },

    // ═══════════════════════════════════════════════════════════════════════
    // PRESCRIPTION DRUGS / AMR MARKERS — chemical half-life & hydrolysis
    // k_chem: 0.05–0.20 d⁻¹   α: clamped to 0.10–0.18
    // ═══════════════════════════════════════════════════════════════════════

    AnalyteProfile {
        name: "Amoxicillin", target_marker: "beta-lactam-LC-MSMS", method: "LC-MS/MS",
        category:     AnalyteCategory::PharmaceuticalAMR,
        baseline_log: 1.2,  noise_std: 0.18,
        decay_rate_k: 0.15, z_threshold: 3.0,
    },
    AnalyteProfile {
        name: "Ciprofloxacin", target_marker: "CIP-LC-MSMS", method: "LC-MS/MS",
        category:     AnalyteCategory::PharmaceuticalAMR,
        baseline_log: 1.0,  noise_std: 0.16,
        decay_rate_k: 0.08, z_threshold: 3.0,
    },
    AnalyteProfile {
        // AMR gene shedding — proxy for community β-lactam resistance burden
        name: "AMR_blaTEM", target_marker: "blaTEM", method: "ddPCR",
        category:     AnalyteCategory::PharmaceuticalAMR,
        baseline_log: 4.5,  noise_std: 0.20,
        decay_rate_k: 0.20, z_threshold: 3.0,
    },
    AnalyteProfile {
        name: "Fluoxetine", target_marker: "FLX-LC-MSMS", method: "LC-MS/MS",
        category:     AnalyteCategory::PharmaceuticalAMR,
        baseline_log: 0.8,  noise_std: 0.15,
        decay_rate_k: 0.05, z_threshold: 3.0,
    },

    // ═══════════════════════════════════════════════════════════════════════
    // ILLICIT SUBSTANCES — metabolic product stability kinetics
    // k_metabolite: 0.20–0.28 d⁻¹   α: 0.18–0.24
    // ═══════════════════════════════════════════════════════════════════════

    AnalyteProfile {
        // Detected as urinary metabolite norfentanyl (more stable than parent)
        name: "Fentanyl", target_marker: "norfentanyl-LC-MSMS", method: "LC-MS/MS",
        category:     AnalyteCategory::IllicitSubstance,
        baseline_log: -0.5, noise_std: 0.30,
        decay_rate_k: 0.25, z_threshold: 2.8,
    },
    AnalyteProfile {
        // Primary urinary metabolite benzoylecgonine is highly stable
        name: "Cocaine", target_marker: "benzoylecgonine-LC-MSMS", method: "LC-MS/MS",
        category:     AnalyteCategory::IllicitSubstance,
        baseline_log:  0.2, noise_std: 0.25,
        decay_rate_k: 0.28, z_threshold: 2.8,
    },
    AnalyteProfile {
        name: "Methamphetamine", target_marker: "MAMP-LC-MSMS", method: "LC-MS/MS",
        category:     AnalyteCategory::IllicitSubstance,
        baseline_log: -0.2, noise_std: 0.28,
        decay_rate_k: 0.22, z_threshold: 2.8,
    },
    AnalyteProfile {
        // Detected as primary metabolite 3,4-HHMA for stability
        name: "MDMA", target_marker: "HHMA-LC-MSMS", method: "LC-MS/MS",
        category:     AnalyteCategory::IllicitSubstance,
        baseline_log: -0.4, noise_std: 0.22,
        decay_rate_k: 0.20, z_threshold: 2.8,
    },

    // ═══════════════════════════════════════════════════════════════════════
    // INDUSTRIAL TOXICANTS — near-zero or adsorption-dominated decay
    // k ≈ 0.002–0.04 d⁻¹   α: floored to 0.10
    // ═══════════════════════════════════════════════════════════════════════

    AnalyteProfile {
        // Lead — conservative tracer; adsorbs to biofilm / settleable solids
        name: "Lead_Pb", target_marker: "Pb-ICP-MS", method: "ICP-MS",
        category:     AnalyteCategory::IndustrialToxicant,
        baseline_log:  0.5, noise_std: 0.12,
        decay_rate_k: 0.01, z_threshold: 3.5,
    },
    AnalyteProfile {
        // PFOA — fluorinated surfactant; essentially non-degrading in sewer
        name: "PFAS_PFOA", target_marker: "PFOA-LC-MSMS", method: "LC-MS/MS",
        category:     AnalyteCategory::IndustrialToxicant,
        baseline_log:  0.3, noise_std: 0.10,
        decay_rate_k: 0.005, z_threshold: 3.5,
    },
    AnalyteProfile {
        // TCE — industrial solvent; slow hydrolysis but partial volatilisation
        name: "Trichloroethylene", target_marker: "TCE-GC-MS", method: "GC-MS",
        category:     AnalyteCategory::IndustrialToxicant,
        baseline_log:  0.2, noise_std: 0.15,
        decay_rate_k: 0.04, z_threshold: 3.5,
    },
    AnalyteProfile {
        // Mercury — negligible aqueous decay; settles with particulates
        name: "Mercury_Hg", target_marker: "Hg-CV-AAS", method: "CV-AAS",
        category:     AnalyteCategory::IndustrialToxicant,
        baseline_log: -0.2, noise_std: 0.12,
        decay_rate_k: 0.002, z_threshold: 3.5,
    },
];

/// Return the subset of the catalog belonging to a given category.
pub fn analytes_by_category(cat: &AnalyteCategory) -> impl Iterator<Item = &'static AnalyteProfile> + '_ {
    ANALYTE_CATALOG.iter().filter(move |a| &a.category == cat)
}

/// Look up a profile by analyte name (case-sensitive).
pub fn analyte_by_name(name: &str) -> Option<&'static AnalyteProfile> {
    ANALYTE_CATALOG.iter().find(|a| a.name == name)
}
