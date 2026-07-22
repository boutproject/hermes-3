import pathlib
import re
import subprocess

requirements_path = pathlib.Path("ci/requirements.txt")
content = requirements_path.read_text()

targets = [
    ("boutdata", "https://github.com/boutproject/boutdata.git", "master"),
    ("xbout", "https://github.com/boutproject/xbout.git", "master"),
    ("xhermes", "https://github.com/boutproject/xhermes.git", "main"),
]

updated = content
changed = False

for package, repo_url, branch in targets:
    latest_sha = subprocess.check_output(
        ["git", "ls-remote", repo_url, f"refs/heads/{branch}"],
        text=True,
    ).split()[0]
    pattern = rf"^{package} @ git\+{re.escape(repo_url)}@[0-9a-f]{{40}}$"
    replacement = f"{package} @ git+{repo_url}@{latest_sha}"
    new_updated, count = re.subn(pattern, replacement, updated, flags=re.MULTILINE)
    if count != 1:
        raise SystemExit(f"Expected exactly one match for {package}, found {count}")
    changed = changed or new_updated != updated
    updated = new_updated

if changed:
    requirements_path.write_text(updated.rstrip("\n") + "\n")
