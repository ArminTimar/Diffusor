"""Stage the documentation for the GitHub Pages site.

Copies the README, the manual, the changelog and the validation report into one
folder that Jekyll then builds (see .github/workflows/pages.yml), together with
the site's layout, stylesheet and script from this folder. Every staged page
gets a title and a description; the README becomes the home page, whose title,
badges and first paragraph move into the header section of the layout. The
navigation is written to _data/nav.yml. Relative links are rewritten so that
they keep working on the site: links to another staged page point to its .html
file, images next to a staged page are copied with it, and links to anything
else (code, data, the licence) point to the file on GitHub.

Usage: python .github/pages/build_site.py <output folder>
"""

import posixpath
import re
import shutil
import sys
from pathlib import Path, PurePosixPath

ROOT = Path(__file__).resolve().parents[2]
HERE = Path(__file__).resolve().parent
REPO_URL = "https://github.com/ArminTimar/Diffusor"

# source path in the repository -> path on the site
PAGES = {
    "README.md": "index.md",
    "manual/USER_GUIDE.md": "manual/USER_GUIDE.md",
    "manual/METHODS.md": "manual/METHODS.md",
    "manual/VERIFICATION.md": "manual/VERIFICATION.md",
    "CHANGELOG.md": "CHANGELOG.md",
    "examples/validation/report.md": "examples/validation/report.md",
}

# (short label for the top bar, long label for the footer, staged page)
NAV = [
    ("Home", "Overview", "index.md"),
    ("Guide", "User guide", "manual/USER_GUIDE.md"),
    ("Methods", "Methods", "manual/METHODS.md"),
    ("Verification", "Verification", "manual/VERIFICATION.md"),
    ("Validation", "Validation report", "examples/validation/report.md"),
    ("Changelog", "Changelog", "CHANGELOG.md"),
]

IMAGE_SUFFIXES = {".png", ".jpg", ".jpeg", ".gif", ".svg", ".webp"}
LINK = re.compile(r"(!?\[[^\]]*\]\()([^)\s]+)(\))")


def site_html(site_md):
    """Site path of the built page for a staged .md path."""
    return site_md[:-3] + ".html"


def _relpath(target, start_dir):
    """Relative POSIX path from the folder start_dir to target."""
    return posixpath.relpath(target, start_dir or ".")


def rewrite_links(text, source, site_path, out):
    source_dir = PurePosixPath(source).parent
    site_dir = PurePosixPath(site_path).parent.as_posix()
    if site_dir == ".":
        site_dir = ""

    def fix(match):
        opening, target, closing = match.groups()
        if re.match(r"^[a-z]+:", target) or target.startswith("#"):
            return match.group(0)
        path, _, anchor = target.partition("#")
        anchor = "#" + anchor if anchor else ""
        repo_path = posixpath.normpath((source_dir / path).as_posix())
        if repo_path in PAGES:
            new = _relpath(site_html(PAGES[repo_path]), site_dir) + anchor
        elif (Path(repo_path).suffix.lower() in IMAGE_SUFFIXES
              and (ROOT / repo_path).is_file()):
            # copy the image next to where the page puts it on the site
            dest = out / posixpath.normpath(
                posixpath.join(site_dir, _relpath(repo_path, source_dir.as_posix())))
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(ROOT / repo_path, dest)
            return match.group(0)
        else:
            kind = "tree" if (ROOT / repo_path).is_dir() else "blob"
            new = f"{REPO_URL}/{kind}/main/{repo_path}{anchor}"
        return opening + new + closing

    return LINK.sub(fix, text)


def plain(text):
    text = re.sub(r"!?\[([^\]]*)\]\([^)]*\)", r"\1", text)
    text = re.sub(r"[*_`]", "", text)
    return " ".join(text.split())


def title_and_description(text):
    lines = text.splitlines()
    title = next((l[2:].strip() for l in lines if l.startswith("# ")), "Diffusor")
    paragraph, started = [], False
    for line in lines:
        stripped = line.strip()
        if not started:
            if stripped.startswith("# "):
                started = True
            continue
        if not stripped:
            if paragraph:
                break
            continue
        if stripped.startswith(("[![", "#", ">", "|", "-", "*", "```")):
            if paragraph:
                break
            continue
        paragraph.append(stripped)
    description = plain(" ".join(paragraph))
    if len(description) > 160:
        description = description[:157].rsplit(" ", 1)[0] + "..."
    return title, description


def split_head(text):
    """Split the README into its first paragraph and the text after it.

    The title line and the badge lines are dropped (the header section of the
    home page shows the name, the DOI and the version). The first two
    sentences of the first paragraph become the header's lead; the rest of
    that paragraph stays as the first paragraph of the body.
    """
    lines = text.splitlines()
    i = next(n for n, l in enumerate(lines) if l.startswith("# ")) + 1
    while i < len(lines) and (not lines[i].strip() or lines[i].lstrip().startswith("[![")):
        i += 1
    paragraph = []
    while i < len(lines) and lines[i].strip():
        paragraph.append(lines[i].strip())
        i += 1
    sentences = re.split(r"(?<=\.)\s+(?=[A-Z])", " ".join(paragraph))
    lead = plain(" ".join(sentences[:2]))
    intro = " ".join(sentences[2:])
    body = "\n".join(lines[i:]).lstrip("\n") + "\n"
    return lead, (intro + "\n\n" + body) if intro else body


def version():
    text = (ROOT / "diffusor/__init__.py").read_text(encoding="utf-8")
    return re.search(r'__version__\s*=\s*"([^"]+)"', text).group(1)


def yaml_string(value):
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def write_nav(out):
    lines = []
    for short, long, target in NAV:
        url = "/" if target == "index.md" else "/" + site_html(target)
        lines.append(f"- title: {yaml_string(short)}\n  long: {yaml_string(long)}\n"
                     f"  url: {yaml_string(url)}\n  home: {str(target == 'index.md').lower()}")
    (out / "_data").mkdir(exist_ok=True)
    (out / "_data/nav.yml").write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def main(out):
    out = Path(out).resolve()
    if out.exists():
        shutil.rmtree(out)
    # layout, stylesheet, script and _config.yml, without this script
    shutil.copytree(HERE, out, ignore=shutil.ignore_patterns("*.py", "__pycache__"))
    # the layout links /favicon.ico on every page; use the program's icon
    shutil.copy2(ROOT / "diffusor/gui/icons/diffusor.ico", out / "favicon.ico")
    write_nav(out)
    for source, site_path in PAGES.items():
        text = (ROOT / source).read_text(encoding="utf-8")
        title, description = title_and_description(text)
        extra = ""
        if source == "README.md":
            title = "Diffusor: diffusion chronometry software"
            lead, text = split_head(text)
            extra = f"hero: true\nlead: {yaml_string(lead)}\nversion: {yaml_string(version())}\n"
        text = rewrite_links(text, source, site_path, out)
        front = (
            "---\n"
            f"title: {yaml_string(title)}\n"
            f"description: {yaml_string(description)}\n"
            f"{extra}"
            "---\n\n"
        )
        dest = out / site_path
        dest.parent.mkdir(parents=True, exist_ok=True)
        dest.write_text(front + text, encoding="utf-8", newline="\n")
        print(f"{source} -> {site_path}")


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "_pages_src")
