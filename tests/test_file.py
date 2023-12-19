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
                self.data_dir / "eqref_OOB.json"
                # / "jetto.eqdsk_out",
            ),
            # volt_seconds_per_radian=True,
            # clockwise_phi=True,
            # to_cocos_index=None,
        )
        print(eqd.cocos)
        a = 1
