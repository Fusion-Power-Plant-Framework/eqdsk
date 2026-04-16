# Using the CLI

This package exposes a command line interface (CLI) that can be used to do a few useful things.

It is available as the `eqdsk` command after pip installing the package.

```bash
pip install eqdsk
eqdsk --help
```

## show

The `show` command can be used to display (print) important information about an eqdsk file.

```bash
eqdsk show path/to/eqdsk/file
```

## plot

The `plot` command can be used to display a plot of the plasma boundary (LCFS) with the magnetic axis centre.

```bash
eqdsk plot path/to/eqdsk/file
```

!!! note
    The `plot` command requires the `matplotlib` package to be installed.

## plot-psi

The `plot-psi` command can be used to display a plot of the poloidal flux (psi) in the eqdsk file.

```bash
eqdsk plot-psi path/to/eqdsk/file
```

!!! note
    The `plot-psi` command requires the `matplotlib` package to be installed.

## convert

The `convert` command can be used to convert an eqdsk file to a different format (.json, .eqdsk or IMAS .nc)
or between COCOS versions.

**Save to .json format and convert between COCOS 1 and 2:**

```bash
eqdsk convert path/to/eqdsk/file -fmt json --from 1 --to 2
```

**Save to .eqdsk format:**

```bash
eqdsk convert path/to/eqdsk/file -fmt eqdsk
```

**Convert between COCOS 1 and 2 and provide the sign of qpsi (if not present in the file):**

```bash
eqdsk convert path/to/eqdsk/file --from 1 --to 2 --qpsi-sign -1
```

## IMAS

The above commands can also be used on an IMAS database. To enable this functionality, the optional `imas` dependencies should be install by running `pip install eqdsk[imas]`.

!!! note
    If converting **from** IMAS to another format, `--from` will be ignored if used because IMAS has a fixed COCOS.
    Likewise, if converting **to** IMAS, then `--to` will be ignored for the same reason.

**Convert an EQDSK file to IMAS:**

In this example, the EQDSK file has COCOS 3. We do not specify a `--to` because IMAS has a fixed COCOS, which is handled within `eqdsk`. The file will be written to the default IMAS format which is a NETCDF (.nc) file.

```bash
eqdsk convert path/to/eqdsk/file -fmt imas --from 3
```

There is also the option to specify the location of the NETCDF file or a URI specifying the HDF5 database.

```bash
eqdsk convert path/to/eqdsk/file -fmt imas --from 3 --imas-uri-out converted/imas/file.nc
```

**Convert a NETCDF IMAS database between I/O data dictionary versions:**

I/0 is split with `:`, in this example `database.nc` has data dictionary version `3.42.0` and `new_database.nc` has data dictionary version `4.0.0`.

```bash
eqdsk convert path/to/imas/database.nc --imas-dd-version 3.42.0:4.0.0 --imas-uri-out path/to/imas/new_database.nc
```

Either filepath could also be replaced with a valid IMAS URI to specify use of the HDF5 backend:

```bash
eqdsk convert imas:hdf5?path=/path/to/hdf5/database --imas-dd-version 3.42.0:4.0.0 --imas-uri-out imas:hdf5?path=/path/to/hdf5/new_database
```
