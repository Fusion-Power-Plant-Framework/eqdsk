from pprint import pprint

from eqdsk.file import EQDSKInterface

EQDSKInterface.DEFAULT_COCOS_CONVENTION = 7
ed = EQDSKInterface.from_file("eudemo.eqdsk")

pprint(ed.cocos_convention)
