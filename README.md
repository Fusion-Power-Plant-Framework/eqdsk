# Eqdsk

[![Hatch project](https://img.shields.io/badge/%F0%9F%A5%9A-Hatch-4051b5.svg)](https://github.com/pypa/hatch)
[![linting - Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)[![Actions status](https://github.com/Fusion-Power-Plant-Framework/eqdsk/actions/workflows/test.yml/badge.svg)](https://github.com//Fusion-Power-Plant-Framework/eqdsk/actions)

An EQDSK reader and writer for GEQDSK (more soon), with COCOS identification and conversion.

There is support for writing an eqdsk to a JSON format (which is now preffered).

We have extended the EQDSK standard to allow for the definition of a CoilSet. ---more on this---

## Setup

This project uses [Hatch](https://hatch.pypa.io/latest/).

Although any python environment manager can be used, we recommend using the default eniornment setup by Hatch.

To start using Hatch, it must be installed and accessible from the command line. See the Hatch [installation](https://hatch.pypa.io/latest/install/) for more.

A simple way to install hatch is to run:

```bash
pip install -g hatch
```

If you can run `hatch -h` then Hatch has been successfully installed.

We recommend setting the `dirs.env` in you hatch config to the following:

```toml
[dirs.env]
virtual = ".hatch"
```

The path to this file can be found by running:

```bash
hatch config find
```

It makes it easier to set the path the the venv in your code editor.

Then run:

```bash
hatch shell
```

This will create the default hatch environment in the project folder. Set your python enviorment path for your editor to `.hatch/eqdsk/bin/python`

## Tests

Run

```bash
hatch run test:tests
```

To run all the tests in a python 3.8 and 3.10 environment
