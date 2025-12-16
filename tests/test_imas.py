# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

from pathlib import Path

import pytest
from imas import DBEntry

from eqdsk.file import IMAS_AVAIL, EQDSKInterface

from ._helpers import compare_dicts

DATA_DIR = Path(__file__).parent / "test_data"

if not IMAS_AVAIL:
    pytest.skip("IMAS is unavailable", allow_module_level=True)


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
def test_imas_write_read(tmp_path, file, cocos, imas_dd_version, coil_comparison):
    """Test an eqdsk file can be read and then written to IMAS"""
    eqdsk = EQDSKInterface.from_file(
        DATA_DIR / file, from_cocos=cocos, qpsi_positive=False
    )

    with DBEntry(tmp_path / "test.nc", "w", dd_version=imas_dd_version) as db:
        eqdsk.write(db, file_format="imas")

    with DBEntry(tmp_path / "test.nc", "r", dd_version=imas_dd_version) as db:
        new_eqdsk = EQDSKInterface.from_imas(db).to_cocos("jetto")

    # IMAS has no defined mechanism for storing coil types
    exclusions = ["coil_types"]

    # Some versions of IMAS cannot store coil geometry
    if not coil_comparison:
        exclusions += ["dzc", "Ic", "dxc", "xc", "zc"]

    assert compare_dicts(
        eqdsk.to_dict(),
        new_eqdsk.to_dict(),
        verbose=True,
        almost_equal=True,
        exclude=exclusions,
    )
