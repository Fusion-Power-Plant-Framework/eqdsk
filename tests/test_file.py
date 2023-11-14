from pathlib import Path

from eqdsk.file import EQDSKInterface


class TestEQDSKInterface:
    """Test the EQDSK interface."""

    @classmethod
    def setup_class(cls):
        cls.data_dir = Path(__file__).parent / "test_data"

    def test_read_strict_geqdsk(self):
        """Read and return the COCOS for the eqdsk."""
        eqd = EQDSKInterface.from_file(
            str(
                self.data_dir
                / "equilibria/EUROfusion_DEMO_2015_equilibria"
                / "Random_Equil_COCOS11.eqdsk",
            ),
            volt_seconds_per_radian=True,
            clockwise_phi=True,
        )

        assert eqd.cocos.index == 11
