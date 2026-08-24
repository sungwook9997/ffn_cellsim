"""Export KB_PRESENTATION.md to a single self-contained HTML file.

No third-party dependency (the ffn_sim env has no markdown lib, and the export
must run anywhere): a small pure-stdlib converter for the markdown subset the
report uses — headings, GFM pipe tables, fenced code, blockquotes, lists, bold/
italic/inline-code, images and links. Every ./figs/*.png is inlined as a base64
data-URI so the result is ONE portable file you can open or present from anywhere.

Run::

    conda activate ffn_sim
    python aleph/outputs/tag_kb/presentation/export_html.py
"""
from __future__ import annotations

import base64
import html
import re
from pathlib import Path

HERE = Path(__file__).parent
SRC = HERE / "KB_PRESENTATION.md"
OUT = HERE / "KB_PRESENTATION.html"

CSS = """
:root { --rag:#6B9AC4; --tag:#3C8D7A; --sys:#5B4B8A; --fail:#C0392B; --pass:#2E8B57; }
* { box-sizing: border-box; }
body { font-family: -apple-system, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
       line-height: 1.6; color: #222; max-width: 980px; margin: 0 auto; padding: 2.5rem 1.5rem 6rem; }
h1 { font-size: 2.0rem; border-bottom: 3px solid var(--sys); padding-bottom: .4rem; }
h2 { font-size: 1.5rem; color: var(--sys); margin-top: 2.6rem; border-top: 1px solid #e3e3ec; padding-top: 1.4rem; }
h3 { font-size: 1.18rem; color: #333; margin-top: 1.6rem; }
p, li { font-size: 1.02rem; }
em { color: #555; }
code { background: #f0f1f6; padding: .12em .4em; border-radius: 4px; font-size: .92em;
       font-family: "SF Mono", Menlo, Consolas, monospace; }
pre { background: #1e2030; color: #e6e8f0; padding: 1rem 1.2rem; border-radius: 8px; overflow-x: auto;
      font-size: .86rem; line-height: 1.45; }
pre code { background: none; padding: 0; color: inherit; }
blockquote { border-left: 4px solid var(--sys); background: #f6f5fb; margin: 1.2rem 0;
             padding: .7rem 1.2rem; color: #333; border-radius: 0 6px 6px 0; }
table { border-collapse: collapse; width: 100%; margin: 1.3rem 0; font-size: .96rem; }
th, td { border: 1px solid #d4d6e2; padding: .55rem .8rem; text-align: left; vertical-align: top; }
th { background: var(--sys); color: #fff; }
tr:nth-child(even) td { background: #f7f8fc; }
img { max-width: 100%; height: auto; display: block; margin: 1.4rem auto;
      border: 1px solid #e3e3ec; border-radius: 8px; box-shadow: 0 2px 10px rgba(0,0,0,.06); }
hr { border: none; border-top: 1px solid #e3e3ec; margin: 2.2rem 0; }
a { color: var(--tag); }
@media print { h2 { page-break-before: always; } body { max-width: none; } }
"""

_CODE_TOKEN = "\x00CODE%d\x00"


def _inline(text: str, codes: list) -> str:
    """Inline markdown -> HTML on an ALREADY html-escaped string."""
    # protect inline code spans (escape their contents, stash, restore last)
    def _stash(m):
        codes.append(html.escape(m.group(1)))
        return _CODE_TOKEN % (len(codes) - 1)
    text = re.sub(r"`([^`]+)`", _stash, text)
    # images: ![alt](src)  -> handled by caller for base64; here links + emphasis
    text = re.sub(r"\[([^\]]+)\]\(([^)]+)\)",
                  lambda m: f'<a href="{m.group(2)}">{m.group(1)}</a>', text)
    text = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", text)
    text = re.sub(r"(?<!\*)\*([^*]+)\*(?!\*)", r"<em>\1</em>", text)
    for i, c in enumerate(codes):
        text = text.replace(_CODE_TOKEN % i, f"<code>{c}</code>")
    return text


def _img_data_uri(src: str) -> str:
    p = (HERE / src.lstrip("./")) if src.startswith("./") else (HERE / src)
    if not p.exists():
        return src  # leave as-is if missing
    b64 = base64.b64encode(p.read_bytes()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def convert(md: str) -> str:
    lines = md.split("\n")
    out: list[str] = []
    i, n = 0, len(lines)
    while i < n:
        line = lines[i]
        stripped = line.strip()

        # fenced code block
        if stripped.startswith("```"):
            i += 1
            buf = []
            while i < n and not lines[i].strip().startswith("```"):
                buf.append(html.escape(lines[i]))
                i += 1
            i += 1  # skip closing fence
            out.append("<pre><code>" + "\n".join(buf) + "</code></pre>")
            continue

        # standalone image line: ![alt](src)
        m = re.match(r"^!\[([^\]]*)\]\(([^)]+)\)\s*$", stripped)
        if m:
            out.append(f'<img alt="{html.escape(m.group(1))}" '
                       f'src="{_img_data_uri(m.group(2))}">')
            i += 1
            continue

        # GFM table: header row + separator
        if stripped.startswith("|") and i + 1 < n and re.match(
                r"^\|[\s:\-|]+\|?\s*$", lines[i + 1].strip()):
            def cells(row):
                return [c.strip() for c in row.strip().strip("|").split("|")]
            header = cells(line)
            i += 2  # skip header + separator
            body = []
            while i < n and lines[i].strip().startswith("|"):
                body.append(cells(lines[i]))
                i += 1
            th = "".join(f"<th>{_inline(html.escape(c), [])}</th>" for c in header)
            trs = []
            for r in body:
                tds = "".join(f"<td>{_inline(html.escape(c), [])}</td>" for c in r)
                trs.append(f"<tr>{tds}</tr>")
            out.append(f"<table><thead><tr>{th}</tr></thead>"
                       f"<tbody>{''.join(trs)}</tbody></table>")
            continue

        # horizontal rule
        if re.match(r"^(-{3,}|\*{3,}|_{3,})$", stripped):
            out.append("<hr>")
            i += 1
            continue

        # heading
        m = re.match(r"^(#{1,6})\s+(.*)$", stripped)
        if m:
            lvl = len(m.group(1))
            out.append(f"<h{lvl}>{_inline(html.escape(m.group(2)), [])}</h{lvl}>")
            i += 1
            continue

        # blockquote
        if stripped.startswith(">"):
            buf = []
            while i < n and lines[i].strip().startswith(">"):
                buf.append(re.sub(r"^\s*>\s?", "", lines[i]))
                i += 1
            out.append("<blockquote>"
                       + _inline(html.escape(" ".join(buf)), []) + "</blockquote>")
            continue

        # list (unordered or ordered)
        if re.match(r"^\s*([-*]|\d+\.)\s+", line):
            ordered = bool(re.match(r"^\s*\d+\.\s+", line))
            tag = "ol" if ordered else "ul"
            items = []
            while i < n and re.match(r"^\s*([-*]|\d+\.)\s+", lines[i]):
                item = re.sub(r"^\s*([-*]|\d+\.)\s+", "", lines[i])
                items.append(f"<li>{_inline(html.escape(item), [])}</li>")
                i += 1
            out.append(f"<{tag}>{''.join(items)}</{tag}>")
            continue

        # blank line
        if not stripped:
            i += 1
            continue

        # paragraph (gather consecutive plain lines)
        buf = []
        while i < n and lines[i].strip() and not re.match(
                r"^(#{1,6}\s|>|\s*([-*]|\d+\.)\s|```|\||!\[)", lines[i].strip()) \
                and not re.match(r"^(-{3,}|\*{3,}|_{3,})$", lines[i].strip()):
            buf.append(lines[i].strip())
            i += 1
        if buf:
            out.append("<p>" + _inline(html.escape(" ".join(buf)), []) + "</p>")
        else:
            i += 1
    return "\n".join(out)


def main() -> None:
    md = SRC.read_text()
    body = convert(md)
    doc = (f"<!doctype html>\n<html lang=\"en\">\n<head>\n<meta charset=\"utf-8\">\n"
           f"<meta name=\"viewport\" content=\"width=device-width, initial-scale=1\">\n"
           f"<title>ffn_cellsim KB — RAG to TAG to Contract-Graph</title>\n"
           f"<style>{CSS}</style>\n</head>\n<body>\n{body}\n</body>\n</html>\n")
    OUT.write_text(doc)
    n_img = doc.count("data:image/png;base64,")
    print(f"wrote {OUT}  ({len(doc) // 1024} KB, {n_img} embedded figures)")


if __name__ == "__main__":
    main()
