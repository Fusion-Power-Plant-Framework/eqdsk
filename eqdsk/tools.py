# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

import operator
import os
import string
from collections.abc import Iterable
from importlib import import_module as imp
from importlib import machinery as imp_mach
from importlib import util as imp_u
from itertools import permutations
from json import JSONEncoder, dumps
from os import listdir
from types import ModuleType
from typing import Any, Tuple, Type, Union

import numpy as np

from eqdsk.log import eqdsk_warn


class NumpyJSONEncoder(JSONEncoder):
    """
    A JSON encoder that can handle numpy arrays.
    """

    def default(self, obj):
        """
        Override the JSONEncoder default object handling behaviour for np.arrays.
        """
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def json_writer(
    data, file=None, return_output=False, *, cls=NumpyJSONEncoder, **kwargs,
):
    """
    Write json in the bluemria style.

    Parameters
    ----------
    data: dict
        dictionary to write to json
    filename: str
        filename to write to
    return_output:bool
        return the json as a string
    cls: JsonEncoder
        json encoder child class
    kwargs: dict
        all further kwargs passed to the json writer

    """
    if file is None and not return_output:
        eqdsk_warn("No json action to take")
        return None

    if "indent" not in kwargs:
        kwargs["indent"] = 4

    the_json = dumps(data, cls=cls, **kwargs)

    if file is not None:
        with open(file, "w") as fh:
            fh.write(the_json)
            fh.write("\n")

    if return_output:
        return the_json


def compare_dicts(
    d1, d2, almost_equal=False, verbose=True, rtol=1e-5, atol=1e-8,
):
    """
    Compares two dictionaries. Will print information about the differences
    between the two to the console. Dictionaries are compared by length, keys,
    and values per common keys

    Parameters
    ----------
    d1: dict
        The reference dictionary
    d2: dict
        The dictionary to be compared with the reference
    almost_equal: bool (default = False)
        Whether or not to use np.isclose and np.allclose for numbers and arrays
    verbose: bool (default = True)
        Whether or not to print to the console
    rtol: float
        The relative tolerance parameter, used if ``almost_eqaul`` is True
    atol: float
        The abosulte tolerance parameter, used if ``almost_eqaul`` is True

    Returns
    -------
    the_same: bool
        Whether or not the dictionaries are the same
    """
    nkey_diff = len(d1) - len(d2)
    k1 = set(d1.keys())
    k2 = set(d2.keys())
    intersect = k1.intersection(k2)
    new_diff = k1 - k2
    old_diff = k2 - k1
    same, different = [], []

    # Define functions to use for comparison in either the array, dict, or
    # numeric cases.
    def dict_eq(value_1, value_2):
        return compare_dicts(
            value_1, value_2, almost_equal, verbose, rtol, atol,
        )

    def array_almost_eq(val1, val2):
        return np.allclose(val1, val2, rtol, atol)

    def num_almost_eq(val1, val2):
        return np.isclose(val1, val2, rtol, atol)

    def array_is_eq(val1, val2):
        return (np.asarray(val1) == np.asarray(val2)).all()

    if almost_equal:
        array_eq = array_almost_eq
        num_eq = num_almost_eq
    else:
        array_eq = array_is_eq
        num_eq = operator.eq

    # Map the comparison functions to the keys based on the type of value in d1.
    comp_map = {
        key: array_eq
        if isinstance(val, (np.ndarray, list))
        else dict_eq
        if isinstance(val, dict)
        else num_eq
        if is_num(val)
        else operator.eq
        for key, val in d1.items()
    }

    # Do the comparison
    for k in intersect:
        v1, v2 = d1[k], d2[k]
        try:
            if comp_map[k](v1, v2):
                same.append(k)
            else:
                different.append(k)
        except ValueError:  # One is an array and the other not
            different.append(k)

    the_same = False
    result = "===========================================================\n"
    if nkey_diff != 0:
        compare = "more" if nkey_diff > 0 else "fewer"
        result += f"d1 has {nkey_diff} {compare} keys than d2" + "\n"
    if new_diff != set():
        result += "d1 has the following keys which d2 does not have:\n"
        new_diff = ["\t" + str(i) for i in new_diff]
        result += "\n".join(new_diff) + "\n"
    if old_diff != set():
        result += "d2 has the following keys which d1 does not have:\n"
        old_diff = ["\t" + str(i) for i in old_diff]
        result += "\n".join(old_diff) + "\n"
    if different:
        result += "the following shared keys have different values:\n"
        different = ["\t" + str(i) for i in different]
        result += "\n".join(different) + "\n"
    if (
        nkey_diff == 0
        and new_diff == set()
        and old_diff == set()
        and different == []
    ):
        the_same = True
    else:
        result += "==========================================================="
        if verbose:
            print(result)
    return the_same


def _get_relpath(folder: str, subfolder: str) -> str:
    path = os.sep.join([folder, subfolder])
    if os.path.isdir(path):
        return path
    else:
        raise ValueError(f"{path} Not a valid folder.")


def get_eqdsk_root() -> str:
    """
    Get the eqdsk root install folder.

    Returns
    -------
    root: str
        The full path to the eqdsk root folder, e.g.:
            '/home/user/code/eqdsk'
    """
    import eqdsk

    path = list(eqdsk.__path__)[0]
    root = os.path.split(path)[0]
    return root


def get_eqdsk_path(path: str = "", subfolder: str = "eqdsk") -> str:
    """
    Get a eqdsk path of a module subfolder. Defaults to root folder.

    Parameters
    ----------
    path: str
        The desired path from which to create a full path
    subfolder: str (default = 'eqdsk')
        The subfolder (from the eqdsk root) in which to create a path
        Defaults to the source code folder, but can be e.g. 'tests', or 'data'

    Returns
    -------
    path: str
        The full path to the desired `path` in the subfolder specified
    """
    root = get_eqdsk_root()
    if "egg" in root:
        return f"/{subfolder}"

    path = path.replace("/", os.sep)
    bpath = _get_relpath(root, subfolder)
    return _get_relpath(bpath, path)


def is_num(thing):
    """
    Determine whether or not the input is a number.

    Parameters
    ----------
    thing: unknown type
        The input which we need to determine is a number or not

    Returns
    -------
    num: bool
        Whether or not the input is a number
    """
    if thing is True or thing is False:
        return False
    if thing is np.nan:
        return False
    try:
        float(thing)
        return True
    except (ValueError, TypeError):
        return False


class NumpyJSONEncoder(JSONEncoder):
    """
    A JSON encoder that can handle numpy arrays.
    """

    def default(self, obj):
        """
        Override the JSONEncoder default object handling behaviour for np.arrays.
        """
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)
