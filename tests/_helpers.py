# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

import operator
from pathlib import Path
from typing import Any

import fortranformat as ff
import numpy as np

from eqdsk.tools import is_num


def get_private_dir():
    """Get private data directory"""
    if (private_data := Path(__file__).parent / "test_data" / "private").is_dir():
        return private_data
    if (
        private_data := Path(__file__).parent.parent.parent / "bluemira-private-data"
    ).is_dir():
        return private_data
    return None


def read_strict_geqdsk(file_path):
    """
    Reads an input EQDSK file in, assuming strict adherence to the
    GEQDSK format. Used to check bluemira outputs can be read by
    external readers.

    Note: The main bluemira GEQDSK reader is more forgiving to
    format variations than this!

    Parameters
    ----------
    file_path: str
        Full path string of the file

    """
    # Create FortranRecordReader objects with the Fortran format
    # edit descriptors to be used to parse the G-EQDSK input.
    f2000 = ff.FortranRecordReader("a48,3i4")
    f2020 = ff.FortranRecordReader("5e16.9")
    f2022 = ff.FortranRecordReader("2i5")
    fCSTM = ff.FortranRecordReader("i5")

    # Define helper function to read in flattened arrays.
    def read_flat_array(fortran_format, array_size):
        """
        Reads in a flat (1D) numpy array from a G-EQDSK file.

        Parameters
        ----------
        fortran_format: ff.FortranRecordReader
            FortranRecordReader object for Fortran format edit descriptor
            to be used to parse the format of each line of the output.
        array_size: int
            Number of elements in array to be read in.

        Returns
        -------
        array: np.array
            1D Numpy array of length array_size populated by elements from
            the GEQDSK file input.
        """
        # Initialise numpy array and read in first line.
        array = np.zeros((array_size,))
        line = fortran_format.read(file.readline())
        # Define a counter to track which column in a line
        # is currently being saved into the array.
        col = 0
        # Populate array. If the column index moves past the
        # end of a line, it is reset to zero and the next line is read.
        for i in range(array_size):
            if col == len(line):
                line = fortran_format.read(file.readline())
                col = 0
            array[i] = line[col]
            col += 1
        return array

    # Open file.
    with open(file_path) as file:
        # Read in data. Variable names are for readability;
        # strict format is as defined at
        # https://w3.pppl.gov/ntcc/TORAY/G_EQDSK.pdf
        _id, _, nx, nz = f2000.read(file.readline())
        _xdim, _zdim, _xcentre, _xgrid1, _zmid = f2020.read(file.readline())
        _xmag, _zmag, _psimag, _psibdry, _bcentre = f2020.read(file.readline())
        _cplasma, _psimag, _, _xmag, _ = f2020.read(file.readline())
        _zmag, _, _psibdry, _, _ = f2020.read(file.readline())
        _fpol = read_flat_array(f2020, nx)
        _pressure = read_flat_array(f2020, nx)
        _ffprime = read_flat_array(f2020, nx)
        _pprime = read_flat_array(f2020, nx)
        _psi = read_flat_array(f2020, nx * nz)
        _qpsi = read_flat_array(f2020, nx)
        nbdry, nlim = f2022.read(file.readline())
        _xbdry_zbdry = read_flat_array(f2020, 2 * nbdry)
        _xlim_zlim = read_flat_array(f2020, 2 * nlim)

        # Read in coil information, as found in the GEQDSK extension
        (ncoil,) = fCSTM.read(file.readline())
        _coil = read_flat_array(f2020, 5 * ncoil)


def compare_dicts(
    d1: dict[str, Any],
    d2: dict[str, Any],
    *,
    almost_equal: bool = False,
    verbose: bool = True,
    rtol: float = 1e-5,
    atol: float = 1e-8,
) -> bool:
    """
    Compares two dictionaries. Will print information about the differences
    between the two to the console. Dictionaries are compared by length, keys,
    and values per common keys

    Parameters
    ----------
    d1:
        The reference dictionary
    d2:
        The dictionary to be compared with the reference
    almost_equal:
        Whether or not to use np.isclose and np.allclose for numbers and arrays
    verbose:
        Whether or not to print to the console
    rtol:
        The relative tolerance parameter, used if ``almost_equal`` is True
    atol:
        The absolute tolerance parameter, used if ``almost_equal`` is True

    Returns
    -------
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
            value_1,
            value_2,
            almost_equal=almost_equal,
            verbose=verbose,
            rtol=rtol,
            atol=atol,
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
        key: (
            array_eq
            if isinstance(val, np.ndarray | list)
            else (
                dict_eq
                if isinstance(val, dict)
                else num_eq
                if is_num(val)
                else operator.eq
            )
        )
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
    if nkey_diff == 0 and new_diff == set() and old_diff == set() and different == []:
        the_same = True
    else:
        result += "==========================================================="
        if verbose:
            print(result)  # noqa: T201
    return the_same
