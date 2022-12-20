# bluemira is an integrated inter-disciplinary design tool for future fusion
# reactors. It incorporates several modules, some of which rely on other
# codes, to carry out a range of typical conceptual fusion reactor design
# activities.
#
# Copyright (C) 2021 M. Coleman, J. Cook, F. Franza, I.A. Maione, S. McIntosh, J. Morris,
#                    D. Short
#
# bluemira is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# bluemira is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with bluemira; if not, see <https://www.gnu.org/licenses/>.

"""
A collection of general helper functions for tests.
"""

import contextlib
import os
from pathlib import Path
from unittest import mock


def combine_text_mock_write_calls(open_mock: mock.MagicMock) -> str:
    """
    Combine the 'write' calls made on a mock created using the
    unittest.mock.mock_open function.

    You may have created a suitable mock object using the following

    .. code-block:: python

        with mock.patch("builtins.open", new_callable=mock.mock_open) as open_mock:
            # write a file

        text_written_to_file = combine_text_mock_write_calls(open_mock)

    """
    write_call_args = open_mock.return_value.write.call_args_list
    if len(write_call_args) == 0:
        return ""
    return "".join([call.args[0] for call in write_call_args])


@contextlib.contextmanager
def file_exists(good_file_path: str, isfile_ref: str):
    """
    Context manager to mock os.path.isfile to return True for a specific
    file path.

    This is useful if the code under test checks for the existence of
    several files, but you want to mock such that only one of those
    files exists.

    .. code-block:: python

        with file_exists("some/file", "module.under.test.os.path.isfile"):
            # do some work pretending "some/file" is a file

    """

    def new_isfile(path: str) -> bool:
        if Path(good_file_path).resolve() == Path(path).resolve():
            return True
        return os.path.exists(path) and not os.path.isdir(path)

    with mock.patch(isfile_ref, new=new_isfile) as is_file_mock:
        yield is_file_mock
