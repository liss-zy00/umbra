"""Smoke import for ``src.main`` (does not execute ``run_pipeline``)."""


def test_main_module_importable():
    import src.main as m

    assert hasattr(m, "run_pipeline")
    assert hasattr(m, "parse_args")
