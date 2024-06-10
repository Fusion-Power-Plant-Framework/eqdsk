# Eqdsk

[![Hatch project](https://img.shields.io/badge/%F0%9F%A5%9A-Hatch-4051b5.svg)](https://github.com/pypa/hatch)
[![linting - Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)
[![Actions status](https://github.com/Fusion-Power-Plant-Framework/eqdsk/actions/workflows/main.yml/badge.svg)](https://github.com//Fusion-Power-Plant-Framework/eqdsk/actions)

An EQDSK reader and writer for GEQDSK (more soon), with COCOS identification and conversion.

There is support for writing an eqdsk to a JSON format (which is now preferred).

We have extended the EQDSK standard to optionally allow for the definition of a CoilSet.

## Setup

We are pip installable therefore for the most recent release:

```bash
pip install eqdsk
```

or for the most recent commit

```bash
pip install git+https://github.com/Fusion-Power-Plant-Framework/eqdsk.git
```

For a developer setup please see [CONTRIBUTING.md](CONTRIBUTING.md#setup-with-hatch)

## Basic Usage

To read in an eqdsk (json or eqdsk) in its raw state:
```python
from eqdsk import EQDSKInterface

EQDSKInterface.from_file('file.json', no_cocos=True)
```
To read in an eqdsk file with a known cocos format and convert it to a given cocos format:
```python
EQDSKInterface.from_file('file.eqdsk', from_cocos=11, to_cocos=17)
```
Alternatively if the direction (clockwise or anticlockwise) and the units of phi (V.s or V.s/rad) are known,
the cocos standard will be calculated for you:
```python
EQDSKInterface.from_file('file.eqdsk', clockwise_phi=True, volt_seconds_per_radian=True)
```
