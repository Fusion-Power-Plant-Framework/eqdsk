from pprint import pprint

from eqdsk.file import EQDSKInterface

ed = EQDSKInterface.from_file("jetto.eqdsk_out")

pprint(ed.cocos_convention)
