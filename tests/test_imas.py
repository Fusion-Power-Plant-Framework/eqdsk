# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

from pathlib import Path

import pytest

from eqdsk.file import IMAS_AVAIL, EQDSKInterface

from ._helpers import compare_dicts

DATA_DIR = Path(__file__).parent / "test_data"

if not IMAS_AVAIL:
    pytest.skip("IMAS is unavailable", allow_module_level=True)
else:
    from imas import DBEntry


@pytest.fixture
def _imas_netcdf_uri(tmp_path):
    return tmp_path / "test.nc"


@pytest.fixture
def _imas_hdf5_uri(tmp_path):
    return f"imas:hdf5?path={(tmp_path / 'imasdb').as_posix()}"


@pytest.mark.parametrize(
    ("file", "cocos"),
    [("jetto.eqdsk_out", "jetto"), ("DN-DEMO_eqref_withCoilNames.json", 3)],
)
@pytest.mark.parametrize(
    ("imas_dd_version", "coil_comparison"),
    [
        (None, True),
        ("3.42.0", False),
        ("4.0.0", False),
        ("4.1.0", True),
    ],
)
@pytest.mark.parametrize("file_format", ["imas", None])
@pytest.mark.parametrize("imas_uri_fixture", ["_imas_netcdf_uri", "_imas_hdf5_uri"])
def test_imas_write_read(
    imas_uri_fixture, file, cocos, imas_dd_version, coil_comparison, file_format, request
):
    """Test an eqdsk file can be read and then written to IMAS"""
    imas_uri = request.getfixturevalue(imas_uri_fixture)

    eqdsk = EQDSKInterface.from_file(
        DATA_DIR / file, from_cocos=cocos, to_cocos=11, qpsi_positive=False
    )

    with DBEntry(imas_uri, "w", dd_version=imas_dd_version) as db:
        eqdsk.write(db, file_format=file_format)

    with DBEntry(imas_uri, "r", dd_version=imas_dd_version) as db:
        new_eqdsk = EQDSKInterface.from_imas(db).to_cocos(11)

    # IMAS has no defined mechanism for storing coil types
    exclusions = ["coil_types"]

    # Some versions of IMAS cannot store coil geometry
    if not coil_comparison:
        exclusions += ["ncoil", "coil_names", "dzc", "Ic", "dxc", "xc", "zc"]

    assert compare_dicts(
        eqdsk.to_dict(),
        new_eqdsk.to_dict(),
        verbose=True,
        almost_equal=True,
        exclude=exclusions,
    )
