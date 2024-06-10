# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

import filecmp
from pathlib import Path

import pytest

from eqdsk.file import EQDSKInterface


def private_files():
    """Get private files"""
    file_path = Path(__file__).parent / "test_data" / "private" / "equilibria"

    return [
        (p, "eqdsk", 17 if "jetto" in p.as_posix() else 11)
        for p in file_path.rglob("*eqdsk*")
    ] + [(p, "json", 11) for p in file_path.rglob("*json")]


class TestEQDSKInterface:
    """Test the EQDSK interface."""

    @classmethod
    def setup_class(cls):
        cls.data_dir = Path(__file__).parent / "test_data"

    def test_read_strict_geqdsk(self, caplog):
        """Read and return the COCOS for the eqdsk."""
        caplog.set_level("INFO")

        eqd_default = EQDSKInterface.from_file(
            self.data_dir / "jetto.eqdsk_out",
        )

        logs = caplog.records
        assert logs[0].levelname == "WARNING"
        assert "(7, 8, 17, 18)" in logs[0].message
        assert all(lg.levelname == "INFO" for lg in logs[1:])
        assert all(
            int(lg.message.split("COCOS")[-1].strip(".")) == ind
            for lg, ind in zip(logs[1:], [7, 11, 11], strict=True)
        )
        caplog.clear()

        eqd_as_cc_2 = EQDSKInterface.from_file(
            self.data_dir / "jetto.eqdsk_out",
            volt_seconds_per_radian=True,
            clockwise_phi=True,
            to_cocos_index=2,
        )

        logs = caplog.records
        assert all(lg.levelname == "INFO" for lg in logs)
        assert all(
            int(lg.message.split("COCOS")[-1].strip(".")) == ind
            for lg, ind in zip(logs, [8, 2, 2], strict=True)
        )
        caplog.clear()

        eqd_no_cc = EQDSKInterface.from_file(
            self.data_dir / "jetto.eqdsk_out",
            no_cocos=True,
        )

        assert eqd_default.cocos.index == EQDSKInterface.DEFAULT_COCOS_INDEX
        assert eqd_as_cc_2.cocos.index == 2
        assert caplog.records == []

        with pytest.raises(ValueError):  # noqa: PT011
            eqd_no_cc.cocos  # noqa: B018

    @pytest.mark.parametrize(
        ("file", "ftype", "ind"),
        [("jetto.eqdsk_out", "eqdsk", 17), ("DN-DEMO_eqref.json", "json", 11)],
    )
    def test_read_write_doesnt_change_file(self, file, ftype, ind, tmp_path):
        eqd_default = EQDSKInterface.from_file(self.data_dir / file, ind)

        eqd_default.write(tmp_path / "test", file_format=ftype)

        filecmp.cmp(self.data_dir / file, tmp_path / f"test.{ftype}", shallow=False)

    @staticmethod
    @pytest.mark.parametrize(("file", "ftype", "ind"), private_files())
    def test_read_write_doesnt_change_file_private(file, ftype, ind, tmp_path):
        path = tmp_path / "private"
        path.mkdir()
        eqd_default = EQDSKInterface.from_file(file, ind)

        eqd_default.write(path / "test", file_format=ftype)

        filecmp.cmp(file, path / f"test.{ftype}", shallow=False)
