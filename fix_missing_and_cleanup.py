#!/usr/bin/env python3
"""
1. Move all missing unique files from exampleFiles_backup into exampleFiles
2. Delete the backup (only duplicates remain)

Run from your Atom-Probe-Toolbox directory:
    cd "/Users/peterfelfer/Dropbox/research/01_research_Erlangen/GIT repositories/Atom Probe Toolbox/Atom-Probe-Toolbox"
    python3 fix_missing_and_cleanup.py
"""
import os, json, shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
EXAMPLE = os.path.join(SCRIPT_DIR, "exampleFiles")
BACKUP = os.path.join(SCRIPT_DIR, "exampleFiles_backup")
MANIFEST = os.path.join(SCRIPT_DIR, "missing_files_manifest.json")

if not os.path.exists(BACKUP):
    print("No exampleFiles_backup found. Nothing to do.")
    exit(0)
if not os.path.exists(MANIFEST):
    print("missing_files_manifest.json not found.")
    exit(1)

with open(MANIFEST) as f:
    manifest = json.load(f)

def fmt(b):
    if b > 1e9: return f"{b/1e9:.1f} GB"
    if b > 1e6: return f"{b/1e6:.1f} MB"
    if b > 1e3: return f"{b/1e3:.1f} KB"
    return f"{b} B"

# Step 1: Move missing files from backup into exampleFiles
print(f"Moving {sum(len(e['files']) for e in manifest)} missing files into exampleFiles...")
moved = 0
moved_size = 0
errors = 0
new_folders = 0

for entry in manifest:
    folder = os.path.join(EXAMPLE, entry['folder'])
    is_new = not os.path.exists(folder)
    os.makedirs(folder, exist_ok=True)
    if is_new:
        new_folders += 1

    for finfo in entry['files']:
        src = os.path.join(BACKUP, finfo['backup_rel'])
        dst = os.path.join(folder, finfo['name'])

        if os.path.exists(dst):
            continue  # already there
        if not os.path.exists(src):
            print(f"  WARNING: source missing: {finfo['backup_rel']}")
            errors += 1
            continue

        try:
            sz = os.path.getsize(src)
            shutil.move(src, dst)
            moved += 1
            moved_size += sz
            if moved % 20 == 0:
                print(f"  Moved {moved} files ({fmt(moved_size)})...")
        except Exception as e:
            print(f"  ERROR moving {finfo['name']}: {e}")
            errors += 1

    # Create/update README.md
    readme = os.path.join(folder, "README.md")
    if not os.path.exists(readme) and entry.get('description'):
        lines = [f"# {entry['folder']}", ""]
        if entry['description']:
            lines.append(f"**Original location:** `exampleFiles_backup/{entry['description']}`")
            lines.append("")
        data_files = [fi for fi in entry['files'] if fi['type'] == 'data']
        range_files = [fi for fi in entry['files'] if fi['type'] == 'range']
        if data_files:
            lines.append("## Data files")
            for fi in data_files:
                lines.append(f"- `{fi['name']}`")
            lines.append("")
        if range_files:
            lines.append("## Range / figure files")
            for fi in range_files:
                lines.append(f"- `{fi['name']}`")
            lines.append("")
        with open(readme, 'w') as f:
            f.write('\n'.join(lines))

print(f"\nDone moving files:")
print(f"  Files moved:       {moved} ({fmt(moved_size)})")
print(f"  New folders:       {new_folders}")
print(f"  Errors:            {errors}")

# Step 2: Also move any remaining files in backup that are symlinked from exampleFiles
print("\nChecking for symlinks pointing to backup...")
backup_real = os.path.realpath(BACKUP)
extra_moved = 0

for root, dirs, files in os.walk(EXAMPLE):
    for f in files:
        fp = os.path.join(root, f)
        if not os.path.islink(fp):
            continue
        try:
            resolved = os.path.realpath(fp)
        except:
            continue
        if resolved.startswith(backup_real + "/"):
            if os.path.exists(resolved):
                try:
                    sz = os.path.getsize(resolved)
                    os.remove(fp)
                    shutil.move(resolved, fp)
                    extra_moved += 1
                    moved_size += sz
                except Exception as e:
                    print(f"  ERROR moving symlinked file {f}: {e}")

if extra_moved:
    print(f"  Moved {extra_moved} additional symlinked files from backup")

# Step 3: Check what's left in backup
remaining = 0
remaining_size = 0
for root, dirs, files in os.walk(BACKUP):
    for f in files:
        try:
            remaining_size += os.path.getsize(os.path.join(root, f))
            remaining += 1
        except:
            pass

print(f"\nRemaining in backup: {remaining} files ({fmt(remaining_size)})")
if remaining > 0:
    print("These are duplicates of files that exist elsewhere in Dropbox.")
    resp = input(f"Delete exampleFiles_backup/ to free {fmt(remaining_size)}? [y/N] ")
    if resp.strip().lower() == 'y':
        shutil.rmtree(BACKUP)
        print("Backup deleted!")
    else:
        print("Backup kept.")
else:
    print("Backup is empty, removing...")
    shutil.rmtree(BACKUP)
    print("Backup deleted!")

# Step 4: Final integrity summary
total_folders = 0
total_data = 0
total_range = 0
total_readme = 0
broken_links = 0
for root, dirs, files in os.walk(EXAMPLE):
    if root == EXAMPLE:
        total_folders = len(dirs)
    for f in files:
        fp = os.path.join(root, f)
        ext = os.path.splitext(f.lower())[1]
        if os.path.islink(fp) and not os.path.exists(fp):
            broken_links += 1
        elif f == "README.md":
            total_readme += 1
        elif ext in {'.pos', '.epos', '.ato'}:
            total_data += 1
        elif ext in {'.rng', '.rrng', '.xrng', '.fig'}:
            total_range += 1

print(f"\n=== FINAL EXAMPLEFILES INTEGRITY ===")
print(f"  Dataset folders:  {total_folders}")
print(f"  Data files:       {total_data}")
print(f"  Range/fig files:  {total_range}")
print(f"  README files:     {total_readme}")
print(f"  Broken symlinks:  {broken_links}")
