from pathlib import Path

from eqdsk.file import EQDSKInterface


class TestEQDSKInterface:
    """Test the EQDSK interface."""

    @classmethod
    def setup_class(cls):
        cls.data_dir = Path(__file__).parent / "test_data"

    def read_strict_geqdsk(self):
        """Read and return the COCOS for the eqdsk."""
        ed = EQDSKInterface.from_file(
            str(self.data_dir / "DN-DEMO_eqref.json"), volt_seconds_per_radian=True, to_cocos_index=4
        )

        print(ed.cocos)
