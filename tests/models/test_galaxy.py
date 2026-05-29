import pandas as pd

from src.models.galaxy import Galaxy
from src.pipeline.constants import DEFAULT_CONTROL_GALAXIES


def test_galaxy_initialization_and_types():
    mock_df = pd.DataFrame({"radius_kpc": [1.0, 2.0], "observed_velocity_kms": [50.0, 60.0]})
    mock_metadata = {"morphology_class": "Dwarf Irregulars (Im)", "hubble_t_type": 10}
    g0 = DEFAULT_CONTROL_GALAXIES[0]

    galaxy = Galaxy(
        name=g0,
        distance_mpc=4.04,
        right_ascension_degrees=193.52187,
        declination_degrees=27.14964,
        redshift=0.001214,
        rotation_dataframe=mock_df,
        metadata=mock_metadata,
    )

    assert galaxy.name == g0
    assert galaxy.distance_mpc == 4.04
    assert isinstance(galaxy.rotation_dataframe, pd.DataFrame)
    assert not galaxy.rotation_dataframe.empty
    assert galaxy.metadata["hubble_t_type"] == 10


def test_galaxy_identity_fields():
    """Dataclass equality with distinct ``DataFrame`` objects is ambiguous — compare explicit fields."""
    df_a = pd.DataFrame({"x": [1]})
    df_b = pd.DataFrame({"x": [1]})
    g1 = Galaxy(name="A", rotation_dataframe=df_a)
    g2 = Galaxy(name="A", rotation_dataframe=df_b)
    g3 = Galaxy(name="B", rotation_dataframe=df_a)
    assert g1.name == g2.name
    assert g1.name != g3.name
    assert g1.rotation_dataframe is df_a
    assert g2.rotation_dataframe is df_b


def test_rotation_dataframe_is_mutable_reference():
    df = pd.DataFrame({"x": [1]})
    gal = Galaxy(name="M", rotation_dataframe=df)
    gal.rotation_dataframe["x"] = 99
    assert df["x"].iloc[0] == 99


def test_default_empty_metadata_and_frame():
    g = Galaxy(name="Empty")
    assert g.rotation_dataframe.shape == (0, 0)
    assert g.metadata == {}
