"""Keep the docs' copies of repository files in step with the originals.

As in virgil:

* the docs home page is the README: ``docs/index.md`` is generated from
  ``README.md``, with links rewritten for the docs;
* the tutorials are executed notebooks: each ``notebooks/<name>.ipynb``
  becomes ``docs/<name>.md``, with its figures in ``docs/generated/``
  (the rendering follows virgil's ``scripts/sync_tutorial_docs.py``).

The design notes stay in ``design/`` and are linked into ``docs/design/``
file by file (Zensical follows file links, not directory links), as are
``PLAN.md`` and ``PROVENANCE.md``.

    python scripts/sync_docs.py          # update
    python scripts/sync_docs.py --check  # fail if anything is out of date
"""

import base64
import json
import os
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
GENERATED = DOCS / "generated"
REPO = "https://github.com/benjaminpope/virgil-vlti/blob/main"
HEADER = (
    "<!-- AUTO-GENERATED FROM README.md by scripts/sync_docs.py. "
    "Edit README.md, not this file. -->\n"
)
# README links that point outside docs/, rewritten for the docs.
TREE = "https://github.com/benjaminpope/virgil-vlti/tree/main"
LINKS = {
    f"]({REPO}/PLAN.md)": "](PLAN.md)",
    f"]({REPO}/PROVENANCE.md)": "](PROVENANCE.md)",
    f"]({TREE}/design)": "](design/gravity_calibration_review.md)",
    f"]({REPO}/design/gravity_calibration_review.md)": "](design/gravity_calibration_review.md)",
    f"]({REPO}/docs/instruments/gravity.md)": "](instruments/gravity.md)",
    f"]({REPO}/docs/instruments/pionier.md)": "](instruments/pionier.md)",
    f"]({REPO}/docs/instruments/matisse.md)": "](instruments/matisse.md)",
}
# Tutorial notebooks, in the order of the docs' navigation.
TUTORIALS = ("calibration", "pair_test", "templates")

_ANSI = re.compile(r"\x1b\[[0-9;]*[a-zA-Z]")


def index_text():
    text = (ROOT / "README.md").read_text()
    for old, new in LINKS.items():
        text = text.replace(old, new)
    return HEADER + text


def _text(value):
    return "".join(value) if isinstance(value, list) else str(value)


def _output_text(output):
    kind = output.get("output_type")
    if kind == "stream":
        return _text(output.get("text", ""))
    if kind in ("execute_result", "display_data"):
        return _text(output.get("data", {}).get("text/plain", ""))
    if kind == "error":
        return _text(output.get("traceback", ""))
    return ""


def render_notebook(name):
    """``(markdown, {image name: png bytes})`` for ``notebooks/<name>.ipynb``."""
    nb = json.loads((ROOT / "notebooks" / f"{name}.ipynb").read_text())
    lines = [f"<!-- AUTO-GENERATED FROM notebooks/{name}.ipynb by scripts/sync_docs.py. -->"]
    images = {}
    for c, cell in enumerate(nb["cells"], start=1):
        source = _text(cell.get("source", "")).rstrip()
        if not source.strip():
            continue
        if cell["cell_type"] == "markdown":
            lines += [source.strip(), ""]
            continue
        lines += ["```python", source, "```", ""]
        for o, output in enumerate(cell.get("outputs", []), start=1):
            png = output.get("data", {}).get("image/png")
            if png is not None:
                image = f"{name}_cell{c:03d}_out{o:02d}.png"
                images[image] = base64.b64decode(_text(png).replace("\n", ""))
                lines += [f"![{name} output {c}.{o}](generated/{image})", ""]
                continue
            text = _ANSI.sub("", _output_text(output)).rstrip()
            # Progress bars and widgets print nothing useful on a static page.
            if text and not text.startswith(("HBox(", "VBox(")):
                lines += ["```text", text, "```", ""]
    return "\n".join(lines).rstrip() + "\n", images


def wanted_links():
    """{link path in docs/: relative target} for every linked file."""
    links = {DOCS / name: Path("..") / name for name in ("PLAN.md", "PROVENANCE.md")}
    for path in sorted((ROOT / "design").rglob("*")):
        if path.is_file():
            rel = path.relative_to(ROOT)
            depth = len(rel.parts)  # docs/design/... sits one level deeper
            links[DOCS / rel] = Path(*[".."] * depth) / rel
    return links


def _sync_file(path, content, check, stale):
    data = content.encode() if isinstance(content, str) else content
    if path.exists() and not path.is_symlink() and path.read_bytes() == data:
        return
    stale.append(path)
    if not check:
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(data)


def main(check=False):
    stale = []
    _sync_file(DOCS / "index.md", index_text(), check, stale)
    wanted_images = set()
    for name in TUTORIALS:
        markdown, images = render_notebook(name)
        _sync_file(DOCS / f"{name}.md", markdown, check, stale)
        for image, png in images.items():
            wanted_images.add(image)
            _sync_file(GENERATED / image, png, check, stale)
    for old in sorted(GENERATED.glob("*.png")) if GENERATED.exists() else ():
        if old.name not in wanted_images:
            stale.append(old)
            if not check:
                old.unlink()
    for link, target in wanted_links().items():
        ok = link.is_symlink() and Path(os.readlink(link)) == target
        if not ok:
            stale.append(link)
            if not check:
                link.parent.mkdir(parents=True, exist_ok=True)
                if link.exists() or link.is_symlink():
                    link.unlink()
                link.symlink_to(target)
    if check and stale:
        print("Out of date (run scripts/sync_docs.py):", *map(str, stale), sep="\n  ")
        return 1
    return 0


if __name__ == "__main__":
    sys.exit(main(check="--check" in sys.argv))
