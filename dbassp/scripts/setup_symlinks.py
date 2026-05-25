"""
setup_symlinks.py
-----------------
One-time setup script to create Windows directory junctions that redirect
local project directories to Google Drive storage.

Junctions created:
  <project>/data/output   -->  G:\\My Drive\\Doctorado\\Automatizacion\\auto_sqlite\\data
  <project>/presentations  -->  G:\\My Drive\\Doctorado\\Automatizacion\\auto_sqlite\\presentations

Run once (as Administrator or with Developer Mode enabled):
    python scripts/setup_symlinks.py

Re-running is safe -- existing junctions are skipped.
"""

import os
import sys
import shutil
import subprocess

# -- Configuration ------------------------------------------------------------
GOOGLE_DRIVE_BASE = r"G:\My Drive\Doctorado\Automatizacion\auto_sqlite"

LINKS = [
    {
        "local":  os.path.join(os.path.dirname(__file__), "..", "output_db"),
        "target": os.path.join(GOOGLE_DRIVE_BASE, "data", "output_db"),
        "label":  "output_db",
    },
    {
        "local":  os.path.join(os.path.dirname(__file__), "..", "presentations"),
        "target": os.path.join(GOOGLE_DRIVE_BASE, "presentations"),
        "label":  "presentations/",
    },
]


def normalise(path: str) -> str:
    return os.path.normpath(os.path.abspath(path))


def create_junction(local: str, target: str) -> None:
    """Create a Windows directory junction from `local` to `target`."""
    result = subprocess.run(
        ["cmd", "/c", "mklink", "/J", local, target],
        capture_output=True, text=True
    )
    if result.returncode != 0:
        raise OSError(result.stderr.strip() or f"mklink failed for {local}")


def setup():
    errors = []

    # Verify Google Drive base exists
    if not os.path.isdir(GOOGLE_DRIVE_BASE):
        print(f"[ERROR] Google Drive base not found: {GOOGLE_DRIVE_BASE}")
        print("        Make sure Google Drive is mounted and the path exists.")
        sys.exit(1)

    for entry in LINKS:
        local  = normalise(entry["local"])
        target = normalise(entry["target"])
        label  = entry["label"]

        print(f"\n-- {label} ------------------------------------------")
        print(f"   local  : {local}")
        print(f"   target : {target}")

        # Ensure the Google Drive target directory exists
        os.makedirs(target, exist_ok=True)
        print(f"   [ok] Target directory ensured")

        # Check if path already exists
        if os.path.isdir(local):
            real = os.path.realpath(local)
            if normalise(real) == target:
                print(f"   [skip] Junction already exists and points to correct target")
                continue
            elif real != local:
                # It's a junction/symlink but points elsewhere
                print(f"   [warn] Path exists but points elsewhere: {real}")
                print(f"          Remove it manually then re-run this script.")
                errors.append(label)
                continue
            else:
                # Regular directory -- move contents to Google Drive then create junction
                print(f"   Moving existing directory contents to Google Drive target...")
                for item in os.listdir(local):
                    src = os.path.join(local, item)
                    dst = os.path.join(target, item)
                    if not os.path.exists(dst):
                        shutil.move(src, dst)
                        print(f"     moved: {item}")
                    else:
                        print(f"     skipped (already in target): {item}")
                shutil.rmtree(local)
                print(f"   Removed local directory: {local}")

        # Create the junction
        try:
            create_junction(local, target)
            print(f"   [ok] Junction created")
        except OSError as e:
            print(f"   [ERROR] Could not create junction: {e}")
            print(f"           Try running this script as Administrator.")
            errors.append(label)

    print("\n" + "=" * 60)
    if errors:
        print(f"Completed with errors for: {', '.join(errors)}")
        print("Fix the errors above and re-run the script.")
        sys.exit(1)
    else:
        print("All junctions created successfully!")
        print(f"\nFiles written to data/output/ and presentations/ will now")
        print(f"live in Google Drive at:\n  {GOOGLE_DRIVE_BASE}")


if __name__ == "__main__":
    setup()
