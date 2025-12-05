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
def test_imas_write_read(file, cocos):
    """Test an eqdsk file can be read and then written to IMAS"""
    eqdsk = EQDSKInterface.from_file(
        DATA_DIR / file, from_cocos=cocos, qpsi_positive=False
    )

    with DBEntry("test.nc", "w") as db:
        eqdsk.write(db, file_format="imas")

    with DBEntry("test.nc", "r") as db:
        new_eqdsk = EQDSKInterface.from_imas(db).to_cocos("jetto")

    assert compare_dicts(
        eqdsk.to_dict(),
        new_eqdsk.to_dict(),
        verbose=True,
        almost_equal=True,
        exclude=["coil_types"],
    )
