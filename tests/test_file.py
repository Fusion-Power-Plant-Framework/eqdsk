# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
import json
from copy import deepcopy
from io import StringIO
from pathlib import Path
from unittest import mock

import numpy as np
import pytest

from eqdsk.cocos import COCOS
from eqdsk.errors import NoSingleConventionError
from eqdsk.file import EQDSKInterface
from tests._helpers import (
    compare_dicts,
    get_differences_from_capture,
    get_private_dir,
    read_strict_geqdsk,
)


def private_files() -> list[tuple[Path, str, int]]:
    """Get private files

    Returns
    -------
    :
        The list of available private eqdsks their filetype and cocos format
    """
    if (pdir := get_private_dir()) is None:
        return []

    file_path = pdir / "equilibria"

    def _cocos(pth):
        pth = Path(*pth.parts[-2:]).as_posix()
        if "jetto" in pth or "COCOS11" in pth:
            return 11
        if ("STEP" in pth and "BLUEPRINT" not in pth) or "DEMO" in pth:
            return 7
        if "COCOS02" in pth:
            return 2
        return 3

    return [
        *((p, "eqdsk", _cocos(p)) for p in file_path.rglob("*eqdsk*")),
        *((p, "json", _cocos(p)) for p in file_path.rglob("*json")),
    ]


class TestEQDSKInterface:
    """Test the EQDSK interface."""

    @classmethod
    def setup_class(cls):
        cls.data_dir = Path(__file__).parent / "test_data"

    @pytest.mark.parametrize("cocos", [11, "jetto", "C11"])
    def test_read_default_cocos(self, cocos, caplog):
        """Read and return the COCOS for the eqdsk."""
        caplog.set_level("INFO")

        eqd_default = EQDSKInterface.from_file(
            self.data_dir / "jetto.eqdsk_out", from_cocos=cocos
        )

        logs = caplog.records
        assert logs[0].levelname == "INFO"
        assert "COCOS 11" in logs[0].message
        assert eqd_default.cocos.index == EQDSKInterface.DEFAULT_COCOS

    def test_read_cocos_specified(self, caplog):
        eqd_as_cc_2 = EQDSKInterface.from_file(
            self.data_dir / "jetto.eqdsk_out",
            volt_seconds_per_radian=True,
            clockwise_phi=True,
            to_cocos=2,
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
            ("DN-DEMO_eqref.json", "json", 3),
            ("eqref_OOB.json", "json", 7),
            ("DN-DEMO_eqref_withCoilNames.json", "json", 3),
        ],
    )
    def test_read_write_doesnt_change_file(self, file, ftype, ind, tmp_path, capsys):
        eqd_default = EQDSKInterface.from_file(
            self.data_dir / file,
            from_cocos=ind,
            qpsi_positive=False,
        )
        eqd_default_nc = EQDSKInterface.from_file(self.data_dir / file, no_cocos=True)
        if ind != 11:
            assert not compare_dicts(
                eqd_default.to_dict(), eqd_default_nc.to_dict(), verbose=True
            )
            capture = capsys.readouterr()

            x_diff_k = {"ffprime", "pprime", "psi", "psibdry", "psimag"}
            if ind != 7:
                x_diff_k |= {"qpsi", "fpol"}

            assert not get_differences_from_capture(capture.out) ^ x_diff_k

        if ind < 10:
            assert abs(eqd_default.psimag) == pytest.approx(
                eqd_default_nc.psimag * 2 * np.pi
            )
        else:
            assert eqd_default.psimag == pytest.approx(eqd_default_nc.psimag)

        eqd_default.write(tmp_path / "test", file_format=ftype)

        eqd_test = EQDSKInterface.from_file(tmp_path / f"test.{ftype}", no_cocos=True)
        eqd_test_d = eqd_test.to_dict()
        eqd_def = eqd_default.to_dict()

        eqd_test.identify(as_cocos=11)
        assert eqd_test._cocos is COCOS.C11  # noqa: SLF001

        eqd_def.pop("name")
        eqd_test_d.pop("name")
        assert compare_dicts(eqd_def, eqd_test_d, verbose=True)

    def test_read_wrong_cocos_raises_ValueError(self):
        with pytest.raises(ValueError, match="No convention found"):
            EQDSKInterface.from_file(self.data_dir / "jetto.eqdsk_out", 17)

    @pytest.mark.parametrize(
        ("setup_keys", "end_keys"),
        [
            (
                {"from_cocos": 3, "to_cocos": 11, "qpsi_positive": False},
                {"from_cocos": 11},
            ),
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
        eqdsk.write(fname, file_format="eqdsk", strict_spec=True)

        # Check eqdsk is readable by Fortran readers.
        # This demands stricter adherence to the G-EQDSK
        # format than eqdsk's main reader.
        read_strict_geqdsk(fname)

        eqdsk.write(fname, file_format="eqdsk", strict_spec=False)
        d2 = EQDSKInterface.from_file(fname, **end_keys).to_dict()

        # Write data read in from test file into a new JSON
        # file, with the suffix "_temp"
        jname = fname.with_suffix("").with_suffix(".json")
        eqdsk.write(jname, file_format="json")
        d3 = EQDSKInterface.from_file(jname, **end_keys).to_dict()

        d1.pop("name")
        d2.pop("name")
        d3.pop("name")

        # Compare dictionaries to check data hasn't
        # been changed.
        assert compare_dicts(d1, d2, verbose=True)
        assert compare_dicts(d1, d3, verbose=True)
        assert compare_dicts(d2, d3, verbose=True)

    def test_write_with_wrong_length_qsi_raises_ValueError(self, tmp_path):
        file = self.data_dir / "DN-DEMO_eqref.json"

        eqdsk = EQDSKInterface.from_file(
            file, from_cocos=3, to_cocos=11, qpsi_positive=False
        )
        eqdsk.qpsi = np.ones(2)

        with pytest.raises(ValueError, match="the length"):
            eqdsk.write(
                Path(tmp_path, f"{Path(file).stem}_temp.eqdsk"), file_format="eqdsk"
            )

    @staticmethod
    @pytest.mark.private
    @pytest.mark.parametrize(("file", "ftype", "ind"), private_files())
    def test_read_write_doesnt_change_file_private(file, ftype, ind, tmp_path):
        path = tmp_path / "private"
        path.mkdir(exist_ok=True)

        eqd_default = EQDSKInterface.from_file(file, from_cocos=ind, qpsi_positive=False)

        if ind != 2 and "Random" not in file.name:
            assert eqd_default.comment is None, len(eqd_default.comment)
        else:
            assert eqd_default.comment is not None
            assert "From" in eqd_default.comment

        eqd_default_nc = EQDSKInterface.from_file(file, no_cocos=True)
        if ind != 11:
            assert not compare_dicts(
                eqd_default.to_dict(), eqd_default_nc.to_dict(), verbose=True
            )

        eqd_default.write(tmp_path / "test", file_format=ftype, strict_spec=False)

        eqd_test = EQDSKInterface.from_file(tmp_path / f"test.{ftype}", no_cocos=True)
        eqd_test_d = eqd_test.to_dict()
        eqd_def = eqd_default.to_dict()

        eqd_test.identify(as_cocos=11)
        assert eqd_test._cocos is COCOS.C11  # noqa: SLF001

        eqd_def.pop("name")
        eqd_test_d.pop("name")
        assert compare_dicts(eqd_def, eqd_test_d, verbose=True)

    def test_derived_field_is_calculated_if_not_given(self):
        data_file = Path(self.data_dir, "DN-DEMO_eqref.json")
        with open(data_file) as f:
            eudemo_sof_data = json.load(f)

        mod_sof_data = deepcopy(eudemo_sof_data)
        for field in ["x", "z", "psinorm"]:
            del mod_sof_data[field]

        with mock.patch(
            "pathlib.Path.open", return_value=StringIO(json.dumps(mod_sof_data))
        ):
            eqdsk = EQDSKInterface.from_file(
                "/some/file.json", from_cocos=3, qpsi_positive=False
            )

        np.testing.assert_allclose(eqdsk.x, eudemo_sof_data["x"])
        np.testing.assert_allclose(eqdsk.z, eudemo_sof_data["z"])
        # The calculation used for psinorm has changed since the
        # eudemo_sof_data was created - so we can't compare to that in
        # this case.
        np.testing.assert_allclose(
            eqdsk.psinorm, np.linspace(0, 1, len(eudemo_sof_data["fpol"]))
        )

    def test_read_matches_values_in_file(self):
        eq = EQDSKInterface.from_file(
            Path(self.data_dir, "jetto.eqdsk_out"), from_cocos=11
        )

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
            mock.patch("pathlib.Path.open", return_value=StringIO("")),
        ):
            EQDSKInterface.from_file(Path(self.data_dir, "jetto.eqdsk_out"), 11)

        with (
            pytest.raises(OSError, match="Should be at least"),
            mock.patch("pathlib.Path.open", return_value=StringIO(" ")),
        ):
            EQDSKInterface.from_file(Path(self.data_dir, "jetto.eqdsk_out"), 11)

    def test_fail_identify_from_file(self):
        with pytest.raises(NoSingleConventionError):
            EQDSKInterface.from_file(Path(self.data_dir, "jetto.eqdsk_out"))

    def test_failed_cocos_fmt(self):
        with pytest.raises(ValueError, match="not a known"):
            EQDSKInterface.from_file(
                Path(self.data_dir, "jetto.eqdsk_out"), from_cocos="hi"
            )
        with pytest.raises(ValueError, match="not a known"):
            EQDSKInterface.from_file(
                Path(self.data_dir, "jetto.eqdsk_out"), from_cocos=b"hi"
            )
        with pytest.raises(ValueError, match="Convention number"):
            EQDSKInterface.from_file(
                Path(self.data_dir, "jetto.eqdsk_out"), from_cocos=10
            )

    @pytest.mark.parametrize("ftype", ["eqdsk", "json"])
    @pytest.mark.parametrize("comment", [True, False])
    def test_eqdsk_with_comment(self, ftype, comment, tmp_path):
        data = Path(self.data_dir, "jetto.eqdsk_out").read_text()
        extra = "\n    some comment\n    over many lines\n    and stuff"
        data += extra
        with mock.patch("pathlib.Path.open", return_value=StringIO(data)):
            eqdsk = EQDSKInterface.from_file(
                Path(self.data_dir, "jetto.eqdsk_out"), from_cocos=11
            )

        assert eqdsk.comment == extra[1:]

        path = tmp_path / f"test.{ftype}"
        eqdsk.write(path, file_format=ftype, write_comment=comment)

        if ftype == "json":
            with path.open() as file:
                f_data = json.load(file)
            if comment:
                assert f_data["comment"] == extra[1:]
            else:
                assert f_data.get("comment") is None

        elif ftype == "eqdsk":
            f_data = path.read_text()
            if comment:
                assert f_data.endswith(f"\n\n{extra}\n")
                assert not f_data.endswith(f"\n\n\n{extra}\n")
            else:
                assert not f_data.endswith(extra)
