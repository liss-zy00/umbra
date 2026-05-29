from src.pipeline.per_galaxy.morphology import catalog_row_for_galaxy, include_in_scientific_evidence_report, is_photometric_lsb
from src.providers.sparc import SparcMorphologyProvider


def test_catalog_and_lsb_flags(sample_morphology_json):
    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    row = catalog_row_for_galaxy("SYN001", lookup, prov)
    assert row.get("morphology_class") == "Dwarf Irregulars (Im)"
    synth = {"galaxy_name": "SYN001", "morphology_class": "Dwarf Irregulars (Im)"}
    assert include_in_scientific_evidence_report(synth, lookup, prov) is True
    lsb_row = {"galaxy_name": "SYN001", "morphology_class": "Late-type Spirals (Sc, Scd)"}
    assert is_photometric_lsb(lsb_row, lookup, prov) is True
