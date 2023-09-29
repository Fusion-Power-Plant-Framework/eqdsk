from eqdsk.cocos import COCOSHelper
from eqdsk.file import EQDSKInterface
from pprint import pprint


ed = EQDSKInterface.from_file("jetto.eqdsk_out")
a = COCOSHelper.identify_eqdsk(ed)

pprint(a)
