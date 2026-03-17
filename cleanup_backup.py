#!/usr/bin/env python3
"""
Move unique files from exampleFiles_backup into exampleFiles dataset folders,
then delete the backup.

Run from your Atom-Probe-Toolbox directory:
    cd "/Users/peterfelfer/Dropbox/research/01_research_Erlangen/GIT repositories/Atom Probe Toolbox/Atom-Probe-Toolbox"
    python3 cleanup_backup.py
"""
import os, shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(SCRIPT_DIR, "exampleFiles")
BACKUP = os.path.join(SCRIPT_DIR, "exampleFiles_backup")

if not os.path.exists(BACKUP):
    print("No exampleFiles_backup found. Nothing to do.")
    exit(0)

if not os.path.exists(EXAMPLE):
    print("No exampleFiles found. Run reorganize_examples.py first.")
    exit(1)

def fmt(b):
    if b > 1e9: return f"{b/1e9:.1f} GB"
    if b > 1e6: return f"{b/1e6:.1f} MB"
    if b > 1e3: return f"{b/1e3:.1f} KB"
    return f"{b} B"

# Step 1: Find all symlinks in exampleFiles that point into exampleFiles_backup
print("Scanning exampleFiles for symlinks pointing to backup...")
backup_real = os.path.realpath(BACKUP)
moved = 0
moved_size = 0
errors = 0

for root, dirs, files in os.walk(EXAMPLE):
    for f in files:
        fp = os.path.join(root, f)
        if not os.path.islink(fp):
            continue

        # Resolve the symlink to its real path
        try:
            resolved = os.path.realpath(fp)
        except:
            continue

        # Check if it points into the backup
        if not resolved.startswith(backup_real + "/"):
            continue

        # This symlink points to a file in the backup.
        # Move the actual file here, replacing the symlink.
        if os.path.exists(resolved):
            try:
                sz = os.path.getsize(resolved)
                os.remove(fp)  # remove the symlink
                shutil.move(resolved, fp)  # move real file in
                moved += 1
                moved_size += sz
                if moved % 20 == 0:
                    print(f"  Moved {moved} files ({fmt(moved_size)})...")
            except Exception as e:
                print(f"  ERROR moving {f}: {e}")
                errors += 1
        else:
            # Broken link to backup - just remove it
            try:
                os.remove(fp)
            except:
                pass

print(f"\nMoved {moved} unique files into exampleFiles ({fmt(moved_size)})")
if errors:
    print(f"Errors: {errors}")

# Step 2: Check what's left in the backup
remaining = 0
remaining_size = 0
for root, dirs, files in os.walk(BACKUP):
    for f in files:
        fp = os.path.join(root, f)
        try:
            remaining_size += os.path.getsize(fp)
            remaining += 1
        except:
            pass

print(f"\nRemaining in backup: {remaining} files ({fmt(remaining_size)})")

if remaining > 0:
    print("These are duplicates that exist elsewhere in Dropbox.")
    resp = input(f"Delete exampleFiles_backup/ ({fmt(remaining_size)})? [y/N] ")
    if resp.strip().lower() == 'y':
        shutil.rmtree(BACKUP)
        print("Backup deleted.")
    else:
        print("Backup kept.")
else:
    print("Backup is empty, removing it.")
    shutil.rmtree(BACKUP)
    print("Done.")
