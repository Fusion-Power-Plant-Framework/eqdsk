# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>
#
# SPDX-License-Identifier: LGPL-2.1-or-later

from json import JSONEncoder, dumps

import numpy as np

from eqdsk.log import eqdsk_warn


class NumpyJSONEncoder(JSONEncoder):
    """A JSON encoder that can handle numpy arrays."""

    def default(self, obj):
        """Override the JSONEncoder default object handling behaviour
        for np.arrays.
        """
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return super().default(obj)


def json_writer(
    data,
    file=None,
    return_output=False,
    *,
    cls=NumpyJSONEncoder,
    **kwargs,
):
    """Write json in the bluemria style.

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
    return None


def is_num(thing):
    """Determine whether or not the input is a number.

    Parameters
    ----------
    thing: unknown type
        The input which we need to determine is a number or not

    Returns:
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
