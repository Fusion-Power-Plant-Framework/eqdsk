from pprint import pprint

from eqdsk.file import EQDSKInterface

# to test, can change the default convention like this
# EQDSKInterface.DEFAULT_COCOS_CONVENTION = 7
ed = EQDSKInterface.from_file("test_eqdsks/eudemo.eqdsk")

pprint(ed.cocos_convention)
