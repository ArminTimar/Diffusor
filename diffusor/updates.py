"""Find out whether a newer Diffusor has been released on GitHub.

Nothing here touches the interface or changes any file. ``fetch_latest`` asks
GitHub's public API for the newest published release, ``is_newer`` compares its
tag with the running version, and ``due`` says whether the once-a-day startup
check should run. Only the standard library is used, and a failure of any kind
raises ``UpdateError`` with a sentence that can be shown as it is.
"""
from __future__ import annotations

import json
import re
from dataclasses import dataclass
from datetime import datetime, timedelta, timezone
from typing import Optional, Tuple
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

REPO = "ArminTimar/Diffusor"
API_URL = f"https://api.github.com/repos/{REPO}/releases/latest"
RELEASES_URL = f"https://github.com/{REPO}/releases"
TIMEOUT_SECONDS = 6.0
MAX_BYTES = 1_000_000           # a release description is a few kilobytes
CHECK_EVERY = timedelta(hours=20)


class UpdateError(Exception):
    """The check could not be completed; the message says why in plain words."""


@dataclass(frozen=True)
class ReleaseInfo:
    version: str        # "0.2.0": the tag without its leading v
    tag: str            # "v0.2.0"
    url: str            # the release page, always on github.com/<REPO>
    notes: str          # the release description, may be empty


def parse_version(text: str) -> Optional[Tuple[int, ...]]:
    """``"v0.2.0"`` -> ``(0, 2, 0)``; None when the text is not a plain version.

    A pre-release suffix (``0.2.0rc1``) makes the version unusable on purpose:
    the update notice is for finished releases.
    """
    m = re.fullmatch(r"[vV]?(\d+(?:\.\d+)*)", (text or "").strip())
    return tuple(int(p) for p in m.group(1).split(".")) if m else None


def is_newer(latest: str, current: str) -> bool:
    """True when ``latest`` is a later version than ``current``. ``1.0`` equals ``1.0.0``."""
    a, b = parse_version(latest), parse_version(current)
    if a is None or b is None:
        return False
    n = max(len(a), len(b))
    return a + (0,) * (n - len(a)) > b + (0,) * (n - len(b))


def release_from_json(data: dict) -> ReleaseInfo:
    """Read the fields Diffusor needs from GitHub's release record."""
    tag = str(data.get("tag_name") or "").strip()
    if parse_version(tag) is None:
        raise UpdateError(f"The latest release is tagged '{tag}', which is not a version "
                          "number, so it cannot be compared with this one.")
    # The link is opened in the user's browser, so accept it only if it points at
    # this project's own release pages.
    url = str(data.get("html_url") or "")
    if not url.startswith(f"https://github.com/{REPO}/releases/"):
        url = f"{RELEASES_URL}/tag/{quote(tag)}"
    return ReleaseInfo(version=tag.lstrip("vV"), tag=tag, url=url,
                       notes=str(data.get("body") or "").strip())


def fetch_latest(current: str, timeout: float = TIMEOUT_SECONDS) -> ReleaseInfo:
    """The newest published release. Raises ``UpdateError`` if it cannot be read."""
    req = Request(API_URL, headers={"Accept": "application/vnd.github+json",
                                    "User-Agent": f"Diffusor/{current}"})
    try:
        with urlopen(req, timeout=timeout) as resp:
            raw = resp.read(MAX_BYTES + 1)
    except HTTPError as exc:
        if exc.code == 404:
            raise UpdateError("No release has been published yet.") from exc
        if exc.code in (403, 429):
            raise UpdateError("GitHub is limiting requests from this network right now. "
                              "Try again later.") from exc
        raise UpdateError(f"GitHub answered with error {exc.code}.") from exc
    except (URLError, TimeoutError, OSError) as exc:
        raise UpdateError("Could not reach GitHub. Check the internet connection.") from exc
    if len(raw) > MAX_BYTES:
        raise UpdateError("GitHub's answer was larger than expected and was ignored.")
    try:
        data = json.loads(raw.decode("utf-8"))
    except (UnicodeDecodeError, ValueError) as exc:
        raise UpdateError("GitHub's answer could not be read.") from exc
    if not isinstance(data, dict):
        raise UpdateError("GitHub's answer could not be read.")
    return release_from_json(data)


def due(last_check_iso: str, now: Optional[datetime] = None) -> bool:
    """True when the startup check has not succeeded in the last ``CHECK_EVERY``."""
    if not last_check_iso:
        return True
    try:
        last = datetime.fromisoformat(last_check_iso)
    except ValueError:
        return True
    if last.tzinfo is None:
        last = last.replace(tzinfo=timezone.utc)
    now = now or datetime.now(timezone.utc)
    return now - last >= CHECK_EVERY or last > now
