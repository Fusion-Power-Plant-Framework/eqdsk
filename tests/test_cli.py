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


def cli_runner(method, args=None, exit_code=0):
    """Click cli runner with test checks"""  # noqa: DOC201
    if args is None:
        args = []

    runner = CliRunner()
    result = runner.invoke(method, args)
    assert result.exit_code == exit_code
    return result


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
        result = cli_runner(cli.cli, ["show", (self.data_dir / filename).as_posix()])
        assert f"COCOS: {cocos}" in result.output

    @pytest.mark.parametrize("filename", ["jetto.eqdsk_out", "DN-DEMO_eqref.json"])
    def test_plot_eqdsk(self, filename):
        pytest.importorskip("matplotlib")
        result = cli_runner(cli.cli, ["plot", (self.data_dir / filename).as_posix()])
        assert not result.output

    @pytest.mark.parametrize("filename", ["jetto.eqdsk_out", "DN-DEMO_eqref.json"])
    def test_plotpsi_eqdsk(self, filename):
        pytest.importorskip("matplotlib")
        result = cli_runner(cli.cli, ["plot-psi", (self.data_dir / filename).as_posix()])
        assert not result.output

    @pytest.mark.parametrize("to", [5, 17])
    @pytest.mark.parametrize(
        ("filename", "fc"), [("jetto.eqdsk_out", 1), ("DN-DEMO_eqref.json", 3)]
    )
    @staticmethod
    def test_convert_eqdsk_with_from_to(filename, fc, to, data_folder):
        result = cli_runner(
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
        assert (data_folder / new_name).is_file()
        assert not result.output

    @pytest.mark.parametrize("filename", ["jetto.eqdsk_out", "DN-DEMO_eqref.json"])
    @staticmethod
    def test_convert_eqdsk_no_cocos(filename, data_folder):
        result = cli_runner(
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
        assert (data_folder / new_name).is_file()
        assert not result.output

    def test_convert_bad_option(self):
        cli_runner(
            cli.cli,
            [
                "convert",
                "-f",
                "1",
                (self.data_dir / "jetto.eqdsk_out").as_posix(),
            ],
            exit_code=2,
        )

        cli_runner(
            cli.cli,
            [
                "convert",
                "-t",
                "11",
                (self.data_dir / "jetto.eqdsk_out").as_posix(),
            ],
            exit_code=2,
        )
