# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

import filecmp
import json
from copy import deepcopy
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

from eqdsk.file import EQDSKInterface
from tests._helpers import compare_dicts, get_private_dir, read_strict_geqdsk


def private_files():
    """Get private files"""
    if (pdir := get_private_dir()) is None:
        return []

    file_path = pdir / "equilibria"

    def _cocos(pth):
        pth = pth.as_posix()
        if "jetto" in pth or "COCOS11" in pth:
            return 11
        return 17

    return [
        *((p, "eqdsk", _cocos(p)) for p in file_path.rglob("*eqdsk*")),
        *((p, "json", _cocos(p)) for p in file_path.rglob("*json")),
    ]


class TestEQDSKInterface:
    """Test the EQDSK interface."""

    @classmethod
    def setup_class(cls):
        cls.data_dir = Path(__file__).parent / "test_data"

    def test_read_default_cocos(self, caplog):
        """Read and return the COCOS for the eqdsk."""
        caplog.set_level("INFO")

        eqd_default = EQDSKInterface.from_file(
            self.data_dir / "jetto.eqdsk_out",
        )

        logs = caplog.records
        assert logs[0].levelname == "WARNING"
        assert "(1, 2, 11, 12)" in logs[0].message
        assert all(lg.levelname == "INFO" for lg in logs[1:])
        assert all(
            int(lg.message.split("COCOS")[-1].strip(".")) == ind
            for lg, ind in zip(logs[1:], [1, 11, 11], strict=True)
        )

        assert eqd_default.cocos.index == EQDSKInterface.DEFAULT_COCOS_INDEX

    def test_read_cocos_specified(self, caplog):
        eqd_as_cc_2 = EQDSKInterface.from_file(
            self.data_dir / "jetto.eqdsk_out",
            volt_seconds_per_radian=True,
            clockwise_phi=True,
            to_cocos_index=2,
        )

        logs = caplog.records
        assert logs[0].levelname == "INFO"
        assert int(logs[0].message.split("COCOS")[-1].strip(".")) == 2
        assert eqd_as_cc_2.cocos.index == 2

    def test_no_cocos_access_raise_ValueError(self):
        with pytest.raises(ValueError, match="not yet been identified"):
            EQDSKInterface.from_file(  # noqa: B018
                self.data_dir / "jetto.eqdsk_out",
                no_cocos=True,
            ).cocos

    @pytest.mark.parametrize(
        ("file", "ftype", "ind"),
        [
            ("jetto.eqdsk_out", "eqdsk", 11),
            ("DN-DEMO_eqref.json", "json", 17),
            ("eqref_OOB.json", "json", 17),
        ],
    )
    def test_read_write_doesnt_change_file(self, file, ftype, ind, tmp_path):
        eqd_default = EQDSKInterface.from_file(self.data_dir / file, ind)

        eqd_default.write(tmp_path / "test", file_format=ftype)

        filecmp.cmp(self.data_dir / file, tmp_path / f"test.{ftype}", shallow=False)

    def test_read_wrong_cocos_raises_ValueError(self):
        with pytest.raises(ValueError, match="No convention found"):
            EQDSKInterface.from_file(self.data_dir / "jetto.eqdsk_out", 17)

    @pytest.mark.parametrize(
        ("setup_keys", "end_keys"),
        [
            ({"from_cocos_index": 17, "to_cocos_index": 11}, {"from_cocos_index": 11}),
            ({"no_cocos": True}, {"no_cocos": True}),
        ],
    )
    def test_read_strict_geqds(self, setup_keys, end_keys, tmp_path):
        file = self.data_dir / "DN-DEMO_eqref.json"
        # Create EQDSK file interface and read data to a dict

        eqdsk = EQDSKInterface.from_file(file, **setup_keys)
        d1 = eqdsk.to_dict()

        # Write data read in from test file into a new EQDSK
        # file, with the suffix "_temp"
        name = Path(file).stem + "_temp"
        fname = Path(tmp_path, f"{name}.eqdsk")
        eqdsk.write(fname, file_format="eqdsk")
        d2 = eqdsk.to_dict()

        # Check eqdsk is readable by Fortran readers.
        # This demands stricter adherence to the G-EQDSK
        # format than eqdsk's main reader.
        read_strict_geqdsk(fname)

        # Write data read in from test file into a new JSON
        # file, with the suffix "_temp"
        jname = fname.with_suffix("").with_suffix(".json")
        eqdsk.write(jname, file_format="json")
        d3 = EQDSKInterface.from_file(jname, **end_keys).to_dict()

        # Compare dictionaries to check data hasn't
        # been changed.
        assert compare_dicts(d1, d2, verbose=True)
        assert compare_dicts(d1, d3, verbose=True)
        assert compare_dicts(d2, d3, verbose=True)

    def test_write_with_wrong_length_qsi_raises_ValueError(self, tmp_path):
        file = self.data_dir / "DN-DEMO_eqref.json"
        # Create EQDSK file interface and read data to a dict

        eqdsk = EQDSKInterface.from_file(file, from_cocos_index=17, to_cocos_index=11)
        eqdsk.qpsi = np.ones(2)

        with pytest.raises(ValueError, match="the length"):
            eqdsk.write(
                Path(tmp_path, f"{Path(file).stem}_temp.eqdsk"), file_format="eqdsk"
            )

    @staticmethod
    @pytest.mark.private()
    @pytest.mark.parametrize(("file", "ftype", "ind"), private_files())
    def test_read_write_doesnt_change_file_private(file, ftype, ind, tmp_path):
        path = tmp_path / "private"
        path.mkdir(exist_ok=True)
        eqd_default = EQDSKInterface.from_file(file, ind)

        eqd_default.write(path / "test", file_format=ftype)

        filecmp.cmp(file, path / f"test.{ftype}", shallow=False)

    def test_derived_field_is_calculated_if_not_given(self):
        data_file = Path(self.data_dir, "DN-DEMO_eqref.json")
        with open(data_file) as f:
            eudemo_sof_data = json.load(f)

        mod_sof_data = deepcopy(eudemo_sof_data)
        for field in ["x", "z", "psinorm"]:
            del mod_sof_data[field]

        with mock.patch(
            "pathlib.Path.open", new=mock.mock_open(read_data=json.dumps(mod_sof_data))
        ):
            eqdsk = EQDSKInterface.from_file("/some/file.json")

        np.testing.assert_allclose(eqdsk.x, eudemo_sof_data["x"])
        np.testing.assert_allclose(eqdsk.z, eudemo_sof_data["z"])
        # The calculation used for psinorm has changed since the
        # eudemo_sof_data was created - so we can't compare to that in
        # this case.
        np.testing.assert_allclose(
            eqdsk.psinorm, np.linspace(0, 1, len(eudemo_sof_data["fpol"]))
        )

    def test_read_matches_values_in_file(self):
        eq = EQDSKInterface.from_file(Path(self.data_dir, "jetto.eqdsk_out"))

        assert eq.nz == 151
        assert eq.nx == 151
        assert eq.xdim == pytest.approx(3.14981545)
        assert eq.ncoil == 0
        assert eq.xc.size == 0
        assert eq.nbdry == 72
        np.testing.assert_allclose(
            eq.xbdry[:3], [0.399993127e01, 0.399150254e01, 0.396906908e01]
        )
        np.testing.assert_allclose(
            eq.zbdry[-3:], [-0.507187454e00, -0.240712636e00, 0.263892047e-01]
        )

    def test_failed_read_eqdsk(self):
        with (
            pytest.raises(OSError, match="Could not read"),
            mock.patch("pathlib.Path.open", new=mock.mock_open(read_data="")),
        ):
            EQDSKInterface.from_file(Path(self.data_dir, "jetto.eqdsk_out"), 11)

        with (
            pytest.raises(OSError, match="Should be at least"),
            mock.patch("pathlib.Path.open", new=mock.mock_open(read_data=" ")),
        ):
            EQDSKInterface.from_file(Path(self.data_dir, "jetto.eqdsk_out"), 11)
