# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

from pathlib import Path

from eqdsk.file import EQDSKInterface


class TestEQDSKInterface:
    """Test the EQDSK interface."""

    @classmethod
    def setup_class(cls):
        cls.data_dir = Path(__file__).parent / "test_data"

    def test_read_strict_geqdsk(self):
        """Read and return the COCOS for the eqdsk."""
        EQDSKInterface.from_file(
            self.data_dir / "jetto.eqdsk_out",
            # volt_seconds_per_radian=True,
            # clockwise_phi=True,
            # to_cocos_index=None,
        )
