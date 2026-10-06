from diffusor.datasets import _resolve_examples_dir


def test_examples_path_uses_installed_data_directory_when_source_tree_is_absent(tmp_path):
    source = tmp_path / "site-packages" / "examples"
    installed = tmp_path / "prefix" / "share" / "diffusor" / "examples"
    validation = installed / "validation"
    profiles = validation / "profiles"
    profiles.mkdir(parents=True)
    (validation / "manifest.json").write_text("[]", encoding="utf8")
    (profiles / "profile.csv").write_text("Distance_um,Fo_mol\n0,80\n", encoding="utf8")

    resolved = _resolve_examples_dir(source_dir=source, data_dir=tmp_path / "prefix")

    assert resolved == installed
    assert (resolved / "validation" / "manifest.json").is_file()
    assert (resolved / "validation" / "profiles" / "profile.csv").is_file()


def test_examples_path_prefers_checkout_data(tmp_path):
    source = tmp_path / "checkout" / "examples"
    source.mkdir(parents=True)

    assert _resolve_examples_dir(source_dir=source, data_dir=tmp_path / "prefix") == source


def test_examples_path_falls_back_to_user_data_directory(tmp_path):
    user_installed = tmp_path / "user-base" / "share" / "diffusor" / "examples"
    user_installed.mkdir(parents=True)

    resolved = _resolve_examples_dir(
        source_dir=tmp_path / "missing-source",
        data_dir=tmp_path / "missing-prefix",
        user_data_dir=tmp_path / "user-base",
    )

    assert resolved == user_installed
