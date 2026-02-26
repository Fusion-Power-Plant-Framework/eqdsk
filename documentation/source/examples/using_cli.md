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

The `convert` command can be used to convert an eqdsk file to a different format (.json or .eqdsk)
or between COCOS versions.

**Save to .json format and convert between COCOS 1 and 2:**

```bash
eqdsk convert path/to/eqdsk/file --fmt json --from 1 --to 2
```

**Save to .eqdsk format:**

```bash
eqdsk convert path/to/eqdsk/file --fmt eqdsk
```

**Convert between COCOS 1 and 2 and provide the sign of qpsi (if not present in the file):**

```bash
eqdsk convert path/to/eqdsk/file --from 1 --to 2 --qpsi-sign -1
```
