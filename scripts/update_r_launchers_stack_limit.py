#!/usr/bin/env python3
"""Apply/check the R startup rule in current shell launchers; never run jobs."""

import argparse
import datetime
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile


OPTION = "--max-ppsize=500000"
RULE = (
    "- Start every R process in current Slurm launchers and shell submission "
    "helpers with `Rscript --max-ppsize=500000` (or the equivalent option for "
    "a direct `R` invocation), including `--list`, `--dry-run`, and inline "
    "`-e` calls. This is a runtime safeguard against pointer-protection stack "
    "overflow, not an analysis parameter. Keep the option exactly once per "
    "invocation; preserve other startup flags. Apply this rule to new/edited "
    "launchers and verify with `python3 update_r_launchers_stack_limit.py "
    "--check` from scripts/. Leave historical directories untouched."
)
# Command starts at the beginning of a line or after a shell command separator
# or opening subshell. Includes calls such as task_map=$(Rscript ...).
CALL = re.compile(r"(?m)(?P<prefix>^[ \t]*|(?<=[(;|&])[ \t]*)Rscript(?=[ \t])")


def patch_launcher(content):
    def replace(match):
        remaining = content[match.end():]
        if re.match(r"[ \t]+--max-ppsize=500000(?=[ \t\r\n]|$)", remaining):
            return match.group(0)
        if re.match(r"[ \t]+--max-ppsize(?:[= \t])", remaining):
            raise ValueError("Different protection-stack setting requires review")
        return match.group(0) + " " + OPTION

    return CALL.sub(replace, content)


def patch_rules(content):
    if RULE in content:
        return content
    heading = "## Slurm conventions\n"
    if heading not in content:
        return content.rstrip() + "\n\n" + heading + "\n" + RULE + "\n"
    return content.replace(heading, heading + "\n" + RULE + "\n", 1)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--dry-run", action="store_true")
    group.add_argument("--apply", action="store_true")
    group.add_argument("--check", action="store_true")
    args = parser.parse_args()
    root = Path(__file__).resolve().parent
    launchers = sorted(list(root.glob("*.slurm")) + list(root.glob("*.sh")))
    changes = {}
    call_count = 0
    # Plan and syntax-check every candidate before modifying any file.
    for path in launchers:
        if path.is_symlink():
            raise ValueError("Refusing to modify symlink: " + str(path))
        content = path.read_text()
        count = len(list(CALL.finditer(content)))
        call_count += count
        updated = patch_launcher(content)
        if updated != content:
            changes[path] = updated
        with tempfile.NamedTemporaryFile(mode="w", suffix=path.suffix) as handle:
            handle.write(updated)
            handle.flush()
            subprocess.run(["bash", "-n", handle.name], check=True)
    rules_path = root / "AGENTS.md"
    if rules_path.is_symlink():
        raise ValueError("Refusing to modify symlink: " + str(rules_path))
    rules = rules_path.read_text() if rules_path.exists() else "# Project Instructions\n"
    updated_rules = patch_rules(rules)
    if updated_rules != rules:
        changes[rules_path] = updated_rules
    print("Checked {} shell files and {} Rscript calls.".format(len(launchers), call_count))
    for path in changes:
        print("Would update: " + path.name)
    if args.apply and changes:
        stamp = datetime.datetime.now().strftime("%Y%m%d_%H%M%S_%f")
        backup = root / ("launcher_stack_backup_" + stamp)
        backup.mkdir()
        for path in changes:
            if path.exists():
                shutil.copy2(str(path), str(backup / path.name))
        for path, content in changes.items():
            path.write_text(content)
        print("Updated {} files; previous versions saved in {}".format(len(changes), backup))
    elif not changes:
        print("PASS: all launcher calls and the project rule are current.")
    elif args.check:
        return 1
    else:
        print("Dry run only; no files changed. Use --apply after review.")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        sys.exit(2)
