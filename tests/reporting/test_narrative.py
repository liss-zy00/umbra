from src.config.thresholds import (
    CALIBRATION_STATUS_OPTIMAL,
    CALIBRATION_STATUS_SUBOPTIMAL,
)
from src.pipeline.constants import DEFAULT_CONTROL_GALAXIES
from src.reporting.narrative import NarrativeGenerator


def test_narrative_contains_key_sections():
    reliable_rows = [
        {
            "galaxy_name": "CamB",
            "morphology_class": "Dwarf Irregulars (Im)",
            "theoretical_gap_coverage_ratio": 0.855,
            "best_scenario_name": "scale_0.100_ion_0.010",
            "scenario_gradient_scale_kpc": 0.1,
            "scenario_ionization_fraction": 0.01,
            "pearson_correlation": 0.99,
            "calibration_status": CALIBRATION_STATUS_OPTIMAL,
        }
    ]
    md = NarrativeGenerator.generate_evidence_markdown(reliable_rows)
    assert "magnetic pressure" in md.lower()
    assert "**CamB**" in md
    assert "85.5%" in md
    assert "Notable systems" in md
    assert "CamB" in md
    # OPTIMAL status must NOT trigger the fallback warning.
    assert "calibração padrão" not in md


def test_narrative_does_not_spotlight_individual_galaxies():
    """The narrative must not contain any per-object spotlight section, even with empty input."""
    md = NarrativeGenerator.generate_evidence_markdown([])
    for name in DEFAULT_CONTROL_GALAXIES:
        assert name not in md
    assert "Proof-of-concept target" not in md


def test_narrative_emits_fallback_warning_when_calibration_suboptimal():
    """Galaxies whose Stokes calibration fell back to the emergency constant must be flagged."""
    reliable_rows = [
        {
            "galaxy_name": "FallbackGal",
            "morphology_class": "Dwarf Irregulars (Im)",
            "theoretical_gap_coverage_ratio": 1.0,
            "best_scenario_name": "scale_0.100_ion_0.010",
            "scenario_gradient_scale_kpc": 0.1,
            "scenario_ionization_fraction": 0.01,
            "pearson_correlation": 0.9,
            "calibration_status": CALIBRATION_STATUS_SUBOPTIMAL,
        }
    ]
    md = NarrativeGenerator.generate_evidence_markdown(reliable_rows)
    assert "Resultados baseados em calibração padrão" in md
    assert "falha na leitura da mediana local" in md
