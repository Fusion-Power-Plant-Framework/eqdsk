# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
from pathlib import Path
from shutil import copytree

import pytest
from click.testing import CliRunner

from eqdsk import cli


@pytest.fixture(scope="module")
def data_folder(tmp_path_factory):
    """Copy data folder to temporary folder for conversion purposes"""  # noqa: DOC201
    data_dir = Path(__file__).parent / "test_data"
    dd = tmp_path_factory.mktemp("data")
    copytree(data_dir, dd / "data")
    return dd / "data"


class TestCli:
    """Test cli"""

    @classmethod
    def setup_class(cls):
        cls.data_dir = Path(__file__).parent / "test_data"

    @pytest.mark.parametrize(
        ("filename", "cocos"),
        [("jetto.eqdsk_out", "1, 2, 11, 12"), ("DN-DEMO_eqref.json", "7, 8, 17, 18")],
    )
    def test_show_eqdsk(self, filename, cocos):
        runner = CliRunner()
        result = runner.invoke(cli.cli, ["show", (self.data_dir / filename).as_posix()])
        assert result.exit_code == 0
        assert f"COCOS: {cocos}" in result.output

    @pytest.mark.parametrize("filename", ["jetto.eqdsk_out", "DN-DEMO_eqref.json"])
    def test_plot_eqdsk(self, filename):
        pytest.importorskip("matplotlib")
        runner = CliRunner()
        result = runner.invoke(cli.cli, ["plot", (self.data_dir / filename).as_posix()])
        assert result.exit_code == 0
        assert not result.output

    @pytest.mark.parametrize("filename", ["jetto.eqdsk_out", "DN-DEMO_eqref.json"])
    def test_plotpsi_eqdsk(self, filename):
        pytest.importorskip("matplotlib")
        runner = CliRunner()
        result = runner.invoke(
            cli.cli, ["plot-psi", (self.data_dir / filename).as_posix()]
        )
        assert result.exit_code == 0
        assert not result.output

    @pytest.mark.xfail(reason="click bug on v8.2/3 #3110 and #824")
    @pytest.mark.parametrize("to", [5, 17])
    @pytest.mark.parametrize(
        ("filename", "fc"), [("jetto.eqdsk_out", 1), ("DN-DEMO_eqref.json", 3)]
    )
    @staticmethod
    def test_convert_eqdsk_with_from_to(filename, fc, to, data_folder):
        runner = CliRunner()
        result = runner.invoke(
            cli.cli,
            [
                "convert",
                "-fmt",
                "json",
                "-f",
                f"{fc}",
                "-t",
                f"{to}",
                "-q",
                "-1",
                (data_folder / filename).as_posix(),
            ],
        )

        new_name = Path(filename).with_stem(f"{Path(filename).stem}_out")
        new_name = new_name.with_suffix(".json")
        assert result.exit_code == 0
        assert (data_folder / new_name).is_file()
        assert not result.output

    @pytest.mark.xfail(reason="click bug on v8.2/3 #3110 and #824")
    @pytest.mark.parametrize("filename", ["jetto.eqdsk_out", "DN-DEMO_eqref.json"])
    @staticmethod
    def test_convert_eqdsk_no_cocos(filename, data_folder):
        runner = CliRunner()
        result = runner.invoke(
            cli.cli,
            [
                "convert",
                "-fmt",
                "eqdsk",
                (data_folder / filename).as_posix(),
            ],
        )

        new_name = Path(filename).with_stem(f"{Path(filename).stem}_out")
        new_name = new_name.with_suffix(
            ".eqdsk" if Path(filename).suffix == "json" else Path(filename).suffix
        )
        assert result.exit_code == 0
        assert (data_folder / new_name).is_file()
        assert not result.output

    def test_convert_bad_option(self):
        runner = CliRunner()

        result = runner.invoke(
            cli.cli,
            [
                "convert",
                "-f",
                "1",
                (self.data_dir / "jetto.eqdsk_out").as_posix(),
            ],
        )
        assert result.exit_code == 2

        result = runner.invoke(
            cli.cli,
            [
                "convert",
                "-t",
                "11",
                (self.data_dir / "jetto.eqdsk_out").as_posix(),
            ],
        )
        assert result.exit_code == 2
