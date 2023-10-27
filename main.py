from pprint import pprint

from eqdsk.cocos import COCOS, Sign, ZeroOne
from eqdsk.file import EQDSKInterface

# to test, can change the default convention like this
ed = EQDSKInterface.from_file(
    "test_eqdsks/cc11.eqdsk", volt_seconds_per_radian=True, to_cocos_index=4
)

pprint(ed.cocos)
