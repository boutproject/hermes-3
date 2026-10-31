"""Prepare resolved physics options for boutdata's expression-aware reader."""

import re
import tempfile
from contextlib import contextmanager
from pathlib import Path


def _strip_run_metadata(text):
    """Omit BOUT's reserved [run] section, including its raw multiline Git diff.

    BoutOptionsFile treats unbalanced parentheses as expression continuations
    and square brackets as section headers. A raw diff can contain both and
    cause subsequent physics sections to be swallowed or incorrectly nested.
    Actual settings headers start at column zero; diff context/addition/removal
    lines have a prefix, so they cannot terminate this metadata section.
    """
    lines, skip = [], False
    for line in text.splitlines(keepends=True):
        header = re.fullmatch(
            r"\[([^\[\]\r\n]+)\][ \t]*(?:[#;].*)?", line.rstrip("\r\n")
        )
        if header:
            section = header[1].strip().lower()
            skip = section == "run" or section.startswith("run:")
        if not skip:
            lines.append(line)
    return "".join(lines)


@contextmanager
def physics_options_file(path):
    """Yield a parser-safe options file, preserving the original file.

    Physics values retain their resolved expressions and command-line overrides.
    Loading needs no writes in the run directory; the temporary copy is removed
    after the options have been read eagerly by xhermes/xbout.
    """
    path = Path(path)
    original = path.read_text()
    physics = _strip_run_metadata(original)
    if physics == original:
        yield path
    else:
        with tempfile.TemporaryDirectory(prefix="hermes-analysis-options-") as folder:
            clean = Path(folder) / path.name
            clean.write_text(physics)
            yield clean
