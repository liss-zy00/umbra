from src.providers.sparc import SparcMorphologyProvider


def test_normalize_and_lookup(sample_morphology_json):
    prov = SparcMorphologyProvider(sample_morphology_json)
    lookup = prov.load_morphology_lookup()
    key = prov.normalize_galaxy_name("SYN001")
    assert key in lookup
    assert lookup[key]["morphology_class"] == "Dwarf Irregulars (Im)"
    assert lookup[key]["disk_surface_brightness_solar_per_pc2"] == 50


def test_classify_hubble_t():
    assert SparcMorphologyProvider.classify_morphology_from_t(10) == "Dwarf Irregulars (Im)"
    assert SparcMorphologyProvider.classify_morphology_from_t(2) == "Early-type Spirals (Sa, Sab)"
