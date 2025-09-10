# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later
"""Eqdsk tools"""

from dataclasses import fields
from json import JSONEncoder, dumps
from pathlib import Path

import numpy as np
import numpy.typing as npt
from pydantic.fields import FieldInfo


class NumpyJSONEncoder(JSONEncoder):
    """A JSON encoder that can handle numpy arrays."""

    def default(self, obj):
        """Override the JSONEncoder default object handling behaviour
        for np.arrays.

        Returns
        -------
        :
            The object in a format json can handle
        """
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def json_writer(
    data: dict,
    file: str | Path | None = None,
    *,
    cls=NumpyJSONEncoder,
    **kwargs,
) -> str:
    """Write json in the bluemria style.

    Parameters
    ----------
    data:
        dictionary to write to json
    file:
        filename to write to
    cls:
        json encoder child class
    kwargs:
        all further kwargs passed to the json writer

    Returns
    -------
    :
        The JSON string
    """
    if "indent" not in kwargs:
        kwargs["indent"] = 4

    the_json = dumps(data, cls=cls, **kwargs)

    if file is not None:
        file = Path(file)
        if file.suffix != ".json":
            file = file.with_suffix(".json")

        with open(file, "w") as fh:
            fh.write(the_json)
            fh.write("\n")

    return the_json


def is_num(thing) -> bool:
    """Determine whether or not the input is a number.

    Parameters
    ----------
    thing:
        The input which we need to determine is a number or not

    Returns
    -------
    :
        Whether or not the input is a number
    """
    if thing in {True, False}:
        return False
    try:
        thing = floatify(thing)
    except (ValueError, TypeError):
        return False
    else:
        return not np.isnan(thing)


def floatify(x: npt.ArrayLike) -> float:
    """Converts the np array or float into a float by returning
    the first element or the element itself.

    Notes
    -----
    This function aims to avoid numpy warnings for float(x) for >0 rank scalars
    it emulates the functionality of float conversion

    Returns
    -------
    :
        The value as a float

    Raises
    ------
    ValueError
        If array like object has more than 1 element
    TypeError
        If object is None
    """
    if x is None:
        raise TypeError("The argument cannot be None")
    return np.asarray(x, dtype=float).item()


def aliases(dcls):
    """Decorator to add aliases to dataclasses

    Requires specification in the dataclass like so

    ``` py
    from pydantic import AliasChoices, Field
    from pydantic.dataclasses import dataclass

    @aliases
    @dataclass
    class Test:
        id: int = Field(alias=AliasChoices("id", "identification"))
    ```

    Parameters
    ----------
    dcls:
        dataclass to modify for aliases

    Returns
    -------
    :
        Modified dataclass
    """
    for field in fields(dcls):
        if isinstance(field.default, FieldInfo):
            for alias in set(field.default.alias.choices).difference([field.name]):
                setattr(
                    dcls,
                    alias,
                    property(
                        lambda self, name=field.name: getattr(self, name),
                        lambda self, value, name=field.name: setattr(self, name, value),
                        doc=f"'{field.name}' alias",
                    ),
                )
    return dcls
