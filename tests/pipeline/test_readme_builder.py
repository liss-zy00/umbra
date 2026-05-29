from src.pipeline.aggregation.readme_builder import build_readme_manual


def test_readme_mentions_magnetic_pressure():
    text = build_readme_manual()
    assert "magnetic_pressure_acceleration" in text
    assert "python src/main.py" in text
