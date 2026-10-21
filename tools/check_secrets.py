"""Block commits and pushes that contain credentials (usernames + passwords, API keys, tokens, key files, screenshots).

Used by git hooks in both repos (fluidpy and aofd); standard library only, so it runs with any Python:
    python tools/check_secrets.py --staged    # pre-commit: scan what is about to be committed
    python tools/check_secrets.py --push      # pre-push: scan every commit being pushed (git passes refs on stdin)
    python tools/check_secrets.py --all       # audit: scan every tracked file as it is now

It finds three kinds of leak and never prints the secret itself, only file, line and kind:
1. **exact values** of your real secrets, read at run time from %USERPROFILE%\\.cdsapirc, _netrc/.netrc, the
   Copernicus Marine credentials file and the Windows Credential Manager vault "aofd" (if `keyring` is installed);
2. **patterns**: a CDS `key:` UUID, a netrc `machine … password …` line, `password = …`, GitHub/AWS/Google tokens,
   private-key blocks (placeholders such as `<your-password>` are ignored);
3. **files that must never be committed**: .cdsapirc, netrc, .env, credential/secret files, key files, and any new
   image or PDF in the repo root (screenshots are where passwords hide; real figures live in sub-folders).
Exit 1 blocks the commit/push. To override a false positive, fix the line; never use --no-verify for a real secret.
"""
from __future__ import annotations

import configparser
import re
import subprocess
import sys
from pathlib import Path

HOME = Path.home()
ZERO = "0" * 40

PLACEHOLDER = re.compile(r"<[^>]*>|\.\.\.|your[-_ ]?(password|token|key|username)|\bset\b|xxxx|\*\*\*\*", re.I)
PATTERNS = {
    "CDS API key": re.compile(r"\bkey:\s*([0-9]{3,}:)?[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}", re.I),
    "netrc login line": re.compile(r"\bmachine\s+\S+\s+login\s+\S+\s+password\s+(\S+)", re.I),
    "password assignment": re.compile(r"""\b(password|passwd|passphrase)\b\s*[:=]\s*["']?([^\s"'<>]{6,})""", re.I),
    "GitHub token": re.compile(r"\b(ghp|gho|ghu|ghs|ghr|github_pat)_[A-Za-z0-9_]{20,}"),
    "AWS access key": re.compile(r"\bAKIA[0-9A-Z]{16}\b"),
    "Google API key": re.compile(r"\bAIza[0-9A-Za-z_\-]{35}\b"),
    "private key block": re.compile(r"-----BEGIN [A-Z ]*PRIVATE KEY-----"),
}
BAD_NAME = re.compile(
    r"(^|/)(\.cdsapirc|_?\.?netrc|\.env(\..+)?|\.copernicusmarine[^/]*|[^/]*(credential|secret)[^/]*|"
    r"[^/]*\.(pem|key|p12|pfx)|id_(rsa|ed25519)[^/]*)$", re.I)
ALLOWED_FILES = {"tools/check_secrets.py", "tools/secrets_vault.py", "tools/setup_credentials.py",
                 "CREDENTIALS_GUIDE.md", "CREDENTIALS_GUIDE.docx"}  # guides/tools about credentials, no values
ROOT_MEDIA = re.compile(r"^[^/]+\.(png|jpe?g|gif|bmp|webp|heic|pdf)$", re.I)
# files that legitimately *describe* credentials with placeholders (still scanned for exact values)
DOC_SUFFIXES = (".md", ".docx")


def known_secret_values() -> list[str]:
    """Collect the real secret strings on this machine (never printed). Short values are skipped as too generic."""
    vals: set[str] = set()
    rc = HOME / ".cdsapirc"
    if rc.exists():
        for ln in rc.read_text(errors="ignore").splitlines():
            if ln.strip().lower().startswith("key:"):
                vals.add(ln.split(":", 1)[1].strip())
    for name in ("_netrc", ".netrc"):
        f = HOME / name
        if f.exists():
            vals.update(m.group(1) for m in PATTERNS["netrc login line"].finditer(f.read_text(errors="ignore")))
    cm = HOME / ".copernicusmarine" / ".copernicusmarine-credentials"
    if cm.exists():
        cp = configparser.ConfigParser()
        cp.read(cm)
        vals.update(v for s in cp.sections() for k, v in cp[s].items() if "pass" in k)
    try:
        import keyring  # optional: the vault written by tools/secrets_vault.py
        for name in ("CDSAPI_KEY", "EARTHDATA_PASSWORD", "COPERNICUSMARINE_PASSWORD"):
            v = keyring.get_password("aofd", name)
            if v:
                vals.add(v)
    except Exception:
        pass
    return [v for v in vals if len(v) >= 8]


def git(*args: str) -> str:
    return subprocess.run(["git", *args], capture_output=True, text=True, encoding="utf-8", errors="replace").stdout


def added_lines_from_diff(diff: str) -> list[tuple[str, int, str]]:
    """(file, line number, text) for every added line in a unified diff."""
    out, path, lineno = [], "", 0
    for ln in diff.splitlines():
        if ln.startswith("+++ "):
            path = ln[6:] if ln.startswith("+++ b/") else ln[4:]
        elif ln.startswith("@@"):
            m = re.search(r"\+(\d+)", ln)
            lineno = int(m.group(1)) if m else 0
        elif ln.startswith("+") and not ln.startswith("+++"):
            out.append((path, lineno, ln[1:]))
            lineno += 1
        elif not ln.startswith("-"):
            lineno += 1
    return out


def scan(lines: list[tuple[str, int, str]], new_files: list[str]) -> list[str]:
    problems: list[str] = []
    secrets = known_secret_values()
    for f in new_files:
        if f in ALLOWED_FILES:
            continue
        if BAD_NAME.search(f):
            problems.append(f"{f}: credential/key file must never be committed")
        elif ROOT_MEDIA.search(f):
            problems.append(f"{f}: image/PDF in the repo root (screenshots can show passwords) - move or delete it")
    for path, n, text in lines:
        if any(s in text for s in secrets):
            problems.append(f"{path}:{n}: contains one of YOUR REAL secrets (token or password)")
            continue
        for kind, rx in PATTERNS.items():
            m = rx.search(text)
            if not m:
                continue
            if kind in ("netrc login line", "password assignment"):
                value = m.groups()[-1]
                if PLACEHOLDER.search(m.group(0)) or (path.endswith(DOC_SUFFIXES) and PLACEHOLDER.search(text)):
                    continue
                if kind == "password assignment" and re.fullmatch(r"[A-Za-z_][\w.]*(\(.*\))?", value):
                    continue  # code like password=args.password or getpass(...)
            problems.append(f"{path}:{n}: looks like a {kind}")
            break
    return problems


def staged() -> tuple[list[tuple[str, int, str]], list[str]]:
    diff = git("diff", "--cached", "-U0", "--no-color", "--text")
    new = [f for f in git("diff", "--cached", "--name-only", "--diff-filter=AR").splitlines() if f]
    return added_lines_from_diff(diff), new


def pushed(stdin: str) -> tuple[list[tuple[str, int, str]], list[str]]:
    lines, new = [], []
    for row in stdin.splitlines():
        parts = row.split()
        if len(parts) != 4 or parts[1] == ZERO:  # deleting a branch: nothing to scan
            continue
        local_sha, remote_sha = parts[1], parts[3]
        rng = [local_sha, "--not", "--remotes"] if remote_sha == ZERO else [f"{remote_sha}..{local_sha}"]
        lines += added_lines_from_diff(git("log", "-p", "-U0", "--no-color", "--text", "--format=", *rng))
        new += [f for f in git("log", "--name-only", "--diff-filter=AR", "--format=", *rng).splitlines() if f]
    return lines, sorted(set(new))


def tracked_now() -> tuple[list[tuple[str, int, str]], list[str]]:
    files = [f for f in git("ls-files").splitlines() if f]
    lines = []
    for f in files:
        p = Path(f)
        if p.suffix.lower() in {".png", ".jpg", ".jpeg", ".gif", ".pdf", ".docx", ".ipynb", ".zip", ".mp4", ".mp3"}:
            continue
        try:
            lines += [(f, i, t) for i, t in enumerate(p.read_text(encoding="utf-8").splitlines(), 1)]
        except (UnicodeDecodeError, OSError):
            continue
    return lines, files


def main(argv: list[str]) -> int:
    mode = argv[0] if argv else "--staged"
    lines, new = {"--staged": staged, "--all": tracked_now}.get(mode, lambda: pushed(sys.stdin.read()))()
    problems = scan(lines, new)
    if problems:
        print("check_secrets: BLOCKED - possible credentials:", file=sys.stderr)
        for p in problems[:40]:
            print("  " + p, file=sys.stderr)
        print("Remove them (git restore --staged <file>, or edit the line) and try again.", file=sys.stderr)
        return 1
    print(f"check_secrets: OK ({len(lines)} lines, {len(new)} files checked)")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
