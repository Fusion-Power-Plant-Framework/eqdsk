# Using the API

This packages exposes the [EQDSKInterface][eqdsk.file.EQDSKInterface] class, which is the main interface for reading and writing eqdsk files.

Below are a few examples of how to use the API.

!!! note
    In the examples, `file_path` can be a string or a `pathlib.Path` object pointing to an eqdsk file.

## Read by specifying COCOS parameters

This is a "regular" use case, where you want to read an eqdsk file and have its COCOS be automatically identified.

There are two required input parameters that must be specified (and know beforehand):

- `clockwise_phi`: Whether the magnetic field is defined in the clockwise or counter-clockwise direction, when viewed from above down.
- `volt_seconds_per_radian`: Whether the flux function is in units of V.s or V.s/2π.

``` py title="main.py"
from eqdsk import EQDSKInterface
ed = EQDSKInterface.from_file(
    file_path,
    clockwise_phi=True,
    volt_seconds_per_radian=True
)
```

This will perform a conversion to the `DEFAULT_COCOS` value (COCOS 11) if no `to_cocos` parameter is provided.

### Perform a specific conversion

``` py title="main.py"
from eqdsk import EQDSKInterface
ed = EQDSKInterface.from_file(
    file_path,
    clockwise_phi=True,
    volt_seconds_per_radian=True,
    to_cocos=13
)
```

Will convert the eqdsk file to COCOS 13.

Some known COCOS numbers from cetain code are prdovided (see [KnownCOCOS][eqdsk.cocos.KnownCOCOS]).

``` py title="main.py"
from eqdsk import EQDSKInterface
ed = EQDSKInterface.from_file(
    file_path,
    clockwise_phi=True,
    volt_seconds_per_radian=True,
    to_cocos="JETTO"
)
```

Will convert the eqdsk file to the known JETTO COCOS (COCOS 11).

### Perform no conversion

``` py title="main.py"
from eqdsk import EQDSKInterface
ed = EQDSKInterface.from_file(
    file_path,
    clockwise_phi=True,
    volt_seconds_per_radian=True,
    to_cocos=None
)
```

### Specifying the sign of qpsi

The COCOS identification procedure uses the sign of qpsi and thus requires there be qpsi data in the eqdsk file. However, qpsi is not required to be stored by the eqdsk format and thus may be None (null) in some files.

If this is case, the `qpsi_sign` parameter must be specified and will set the qpsi data to an array (+/-)1's.

``` py title="main.py"
from eqdsk import EQDSKInterface
ed = EQDSKInterface.from_file(
    file_path,
    clockwise_phi=True,
    volt_seconds_per_radian=True,
    qpsi_positive=False
)
```

## Read by specifying the COCOS number

The `clockwise_phi` and `volt_seconds_per_radian` parameters do not need to be set if you know the COCOS number already.

``` py title="main.py"
from eqdsk import EQDSKInterface
ed = EQDSKInterface.from_file(
    file_path,
    from_cocos=8 # Could also be a known string, like "JETTO"
)
```

This will only work if the `from_cocos` number is one of the possible identifiable versions.

There are 4 possible versions if the `clockwise_phi` and `volt_seconds_per_radian` parameters are not set (two if one of them is).

## Read with no COCOS identification or conversions

You may want to read an eqdsk file and access the raw data, without any kind of COCOS identification or conversion.

``` py title="main.py"
from eqdsk import EQDSKInterface
ed = EQDSKInterface.from_file(file_path, no_cocos=True)
```
