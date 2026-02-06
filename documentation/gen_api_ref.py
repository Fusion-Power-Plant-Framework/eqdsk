# SPDX-FileCopyrightText: 2023-present The Bluemira Developers <https://github.com/Fusion-Power-Plant-Framework/bluemira>  # noqa: INP001
#
# SPDX-License-Identifier: LGPL-2.1-or-later

"""Generate the API reference pages and navigation."""

from pathlib import Path

import mkdocs_gen_files
from mkdocs_gen_files.nav import Nav

nav = Nav()

root = Path(__file__).parent.parent
src = root / "eqdsk"

# With the navigation.sections feature enabled, this isn't necessary.
# nav["Overview",] = "overview.md"  # noqa: ERA001

for path in sorted(src.rglob("*.py")):
    module_path = path.relative_to(root).with_suffix("")
    doc_path = path.relative_to(root).with_suffix(".md")
    full_doc_path = Path("reference", doc_path)

    parts = tuple(module_path.parts)

    is_init = False

    if parts[-1] == "__init__":
        is_init = True
        parts = parts[:-1]
    elif parts[-1] == "__main__" or parts[-1] == "_version" or parts[-1] == "cli":
        continue

    p = ".".join(parts)
    parts_p = (*parts, p) if is_init else (*parts[:-1], p)
    nav[parts_p] = doc_path.as_posix()

    with mkdocs_gen_files.open(full_doc_path, "w") as fd:
        fd.write(f"::: {p}")

    mkdocs_gen_files.set_edit_path(full_doc_path, path.relative_to(root))

with mkdocs_gen_files.open("reference/overview.md", "w") as nav_file:
    nav_file.writelines(nav.build_literate_nav())
