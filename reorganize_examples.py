#!/usr/bin/env python3
"""
Reorganize exampleFiles into symlinked dataset folders.
Run this from your Atom-Probe-Toolbox directory:
    cd "/Users/peterfelfer/Dropbox/research/01_research_Erlangen/GIT repositories/Atom Probe Toolbox/Atom-Probe-Toolbox"
    python3 reorganize_examples.py
"""
import os, json, shutil

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
TOOLBOX = SCRIPT_DIR
EXAMPLE = os.path.join(TOOLBOX, "exampleFiles")
BACKUP = os.path.join(TOOLBOX, "exampleFiles_backup")
MANIFEST = os.path.join(TOOLBOX, "dataset_manifest.json")

with open(MANIFEST) as f:
    manifest = json.load(f)

# Step 1: Handle existing directories
if os.path.exists(EXAMPLE) and not os.path.islink(EXAMPLE):
    # If there's already a backup, the current exampleFiles is the broken one from the VM
    if os.path.exists(BACKUP):
        print(f"Removing broken exampleFiles from VM attempt...")
        shutil.rmtree(EXAMPLE)
    else:
        print(f"Backing up original exampleFiles -> exampleFiles_backup")
        os.rename(EXAMPLE, BACKUP)

if not os.path.exists(BACKUP):
    print("ERROR: exampleFiles_backup not found! Cannot proceed.")
    exit(1)

print(f"Creating new exampleFiles/ with {len(manifest)} dataset folders...")
os.makedirs(EXAMPLE, exist_ok=True)

created = 0
symlinks = 0
broken_links = 0

def fmt_size(b):
    if b > 1e9: return f"{b/1e9:.1f} GB"
    if b > 1e6: return f"{b/1e6:.1f} MB"
    if b > 1e3: return f"{b/1e3:.1f} KB"
    return f"{b} B"

for ds in manifest:
    folder = os.path.join(EXAMPLE, ds['folder'])
    os.makedirs(folder, exist_ok=True)
    
    # Symlink data files
    for df in ds['data_files']:
        link = os.path.join(folder, df['name'])
        target = df['host_path']
        
        # Check if target exists
        if not os.path.exists(target):
            # Try the backup
            # The target might be in exampleFiles which is now exampleFiles_backup
            if '/exampleFiles/' in target:
                alt = target.replace('/exampleFiles/', '/exampleFiles_backup/')
                if os.path.exists(alt):
                    target = alt
        
        if os.path.exists(target):
            rel = os.path.relpath(target, folder)
            if not os.path.exists(link):
                os.symlink(rel, link)
                symlinks += 1
        else:
            broken_links += 1
    
    # Symlink range files
    for rf in ds['range_files']:
        link = os.path.join(folder, rf['name'])
        target = rf['host_path']
        
        if not os.path.exists(target):
            if '/exampleFiles/' in target:
                alt = target.replace('/exampleFiles/', '/exampleFiles_backup/')
                if os.path.exists(alt):
                    target = alt
        
        if os.path.exists(target):
            rel = os.path.relpath(target, folder)
            if not os.path.exists(link):
                os.symlink(rel, link)
                symlinks += 1
        else:
            broken_links += 1
    
    # Create README.md
    readme = os.path.join(folder, "README.md")
    lines = [f"# {ds['folder']}", ""]
    if ds['description']:
        lines.append(f"**Context:** {ds['description']}")
        lines.append("")
    lines.append("## Data files")
    for df in ds['data_files']:
        lines.append(f"- `{df['name']}` ({fmt_size(df['size'])})")
    lines.append("")
    if ds['range_files']:
        lines.append("## Range files")
        for rf in ds['range_files']:
            lines.append(f"- `{rf['name']}`")
        lines.append("")
    lines.append("## Original location")
    for df in ds['data_files']:
        lines.append(f"- `{df['host_path']}`")
    lines.append("")
    with open(readme, 'w') as fh:
        fh.write("\n".join(lines))
    
    created += 1

print(f"\nDone!")
print(f"  Folders created:  {created}")
print(f"  Symlinks created: {symlinks}")
print(f"  Broken (missing): {broken_links}")
print(f"\nYour original files are safe in exampleFiles_backup/")
print(f"Once you confirm everything works, you can delete exampleFiles_backup/ to free ~81 GB.")
print(f"(But keep it if any symlinks point to it!)")
