"""Reader for the Ansys Electronics Desktop project text format (.aedt / .aedtz).

The format is a tree of ``$begin 'Name' ... $end 'Name'`` blocks containing
``key=value`` lines and function-style entries such as ``Objects(2530)`` or
``VariableProp('a', 'UD', '', '60mm')``. After the main ``AnsoftProject``
block the file may carry binary attachments (e.g. Parasolid ``.x_b`` bodies),
which are skipped.
"""
from __future__ import annotations

import io
import os
import re
import zipfile

_BEGIN = re.compile(r"^\$begin '(.*)'$")
_END = re.compile(r"^\$end '(.*)'$")
_KV = re.compile(r"^('(?:[^']*)'|[^=(\[]+?)=(.*)$")
_CALL = re.compile(r"^([A-Za-z_][\w ]*?)\s*([(\[])(.*)$")


class Block:
    __slots__ = ("name", "props", "entries", "children")

    def __init__(self, name: str):
        self.name = name
        self.props: dict[str, str] = {}
        self.entries: list[tuple[str, str]] = []   # function-style lines, in order
        self.children: list[Block] = []

    # ------------------------------------------------------------------ navigation
    def child(self, name: str) -> "Block | None":
        for c in self.children:
            if c.name == name:
                return c
        return None

    def children_named(self, name: str) -> list["Block"]:
        return [c for c in self.children if c.name == name]

    def find(self, *path: str) -> "Block | None":
        b = self
        for p in path:
            b = b.child(p) if b is not None else None
        return b

    def walk(self):
        stack = [self]
        while stack:
            b = stack.pop()
            yield b
            stack.extend(reversed(b.children))

    # ------------------------------------------------------------------ values
    def get(self, key: str, default=None) -> str | None:
        """String value with surrounding quotes removed."""
        raw = self.props.get(key)
        return default if raw is None else unquote(raw)

    def get_raw(self, key: str, default=None):
        return self.props.get(key, default)

    def calls(self, name: str) -> list[list[str]]:
        """Arguments of every function-style entry with this name."""
        return [split_args(args) for n, args in self.entries if n == name]

    def call(self, name: str) -> list[str] | None:
        c = self.calls(name)
        return c[0] if c else None

    def __repr__(self):
        return f"Block({self.name!r}, {len(self.props)} props, {len(self.children)} children)"


def unquote(raw: str) -> str:
    s = raw.strip()
    if len(s) >= 2 and s[0] == s[-1] and s[0] in "'\"":
        s = s[1:-1]
    if len(s) >= 2 and s[0] == s[-1] == '"':
        s = s[1:-1]
    return s


def split_args(s: str) -> list[str]:
    """Split a function/array argument string at top-level commas.

    ``"(1, 'a,b', f(2, 3))"`` -> ``['1', 'a,b', 'f(2, 3)']``; array syntax
    ``"[2: 5, 6]"`` drops the leading count.
    """
    s = s.strip()
    if s and s[0] in "([":
        close = ")" if s[0] == "(" else "]"
        if s.endswith(close):
            s = s[1:-1]
        else:
            s = s[1:]
    if re.match(r"^\d+:", s):
        s = s.split(":", 1)[1]
    out, depth, cur, quote = [], 0, [], None
    for ch in s:
        if quote:
            cur.append(ch)
            if ch == quote:
                quote = None
            continue
        if ch in "'\"":
            quote = ch
            cur.append(ch)
        elif ch in "([":
            depth += 1
            cur.append(ch)
        elif ch in ")]":
            depth -= 1
            cur.append(ch)
        elif ch == "," and depth == 0:
            out.append("".join(cur).strip())
            cur = []
        else:
            cur.append(ch)
    tail = "".join(cur).strip()
    if tail:
        out.append(tail)
    return [unquote(a) for a in out]


def ints(s: str | list) -> list[int]:
    items = split_args(s) if isinstance(s, str) else s
    out = []
    for a in items:
        try:
            out.append(int(float(a)))
        except ValueError:
            pass
    return out


def floats(items) -> list[float]:
    items = split_args(items) if isinstance(items, str) else items
    return [float(a) for a in items]


def _read_text(path: str) -> str:
    if path.lower().endswith(".aedtz"):
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.lower().endswith(".aedt")]
            if not names:
                raise ValueError(f"{path}: archive contains no .aedt file")
            data = z.read(names[0])
    else:
        with open(path, "rb") as fh:
            data = fh.read()
    # Main tree is ASCII; binary attachments follow it, so decode permissively.
    return data.decode("latin-1")


def parse(path_or_text: str, is_text: bool = False) -> Block:
    """Parse an .aedt/.aedtz file (or text) and return the 'AnsoftProject' block."""
    text = path_or_text if is_text else _read_text(path_or_text)
    root = Block("<root>")
    stack = [root]
    lines = io.StringIO(text)
    in_cdata = False
    for line in lines:
        s = line.strip()
        if in_cdata:
            if "$end_cdata$" in s:
                in_cdata = False
            continue
        if not s:
            continue
        if s.startswith("$begin '"):
            m = _BEGIN.match(s)
            if m:
                b = Block(m.group(1))
                stack[-1].children.append(b)
                stack.append(b)
                continue
        if s.startswith("$end '"):
            if _END.match(s) and len(stack) > 1:
                stack.pop()
                if len(stack) == 1 and root.children and root.children[0].name == "AnsoftProject":
                    break  # main tree done; the rest is binary attachments
                continue
        cur = stack[-1]
        m = _KV.match(s)
        if m and not s.startswith("$"):
            key = unquote(m.group(1))
            val = m.group(2)
            if "$begin_cdata$" in val and "$end_cdata$" not in val:
                in_cdata = True
            cur.props[key] = val
            continue
        m = _CALL.match(s)
        if m:
            cur.entries.append((m.group(1).strip(), m.group(2) + m.group(3)))
    proj = root.child("AnsoftProject")
    if proj is None:
        raise ValueError("Not an Ansys Electronics Desktop project (no 'AnsoftProject' block)")
    return proj


def design_blocks(project: Block) -> list[Block]:
    """All HFSS designs in the project, in file order."""
    return [b for b in project.children if b.name == "HFSSModel"]


def design_names(path: str) -> list[str]:
    return [d.get("Name", f"Design{i + 1}") for i, d in enumerate(design_blocks(parse(path)))]


def describe(path: str) -> str:
    p = parse(path)
    names = [d.get("Name") for d in design_blocks(p)]
    return f"{os.path.basename(path)}: {len(names)} HFSS design(s): {', '.join(names)}"
