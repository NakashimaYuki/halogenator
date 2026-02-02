#!/usr/bin/env python3
"""
Pack input data for transfer to Linux server.
Creates a compressed tar archive of nplike_v2 data.
"""

import os
import tarfile
import sys
from pathlib import Path
from datetime import datetime

def main():
    # Get project root
    script_dir = Path(__file__).parent
    project_root = script_dir.parent

    # Source directory
    source_dir = project_root / "data" / "output" / "nplike_v2"

    if not source_dir.exists():
        print(f"ERROR: Source directory not found: {source_dir}")
        sys.exit(1)

    # Output archive
    timestamp = datetime.now().strftime("%Y%m%d")
    archive_name = f"nplike_v2_data_{timestamp}.tar.gz"
    archive_path = project_root / archive_name

    print(f"=" * 60)
    print(f"Halogenator Data Packer")
    print(f"=" * 60)
    print(f"Source: {source_dir}")
    print(f"Output: {archive_path}")
    print()

    # Calculate source size
    total_size = 0
    file_count = 0
    for root, dirs, files in os.walk(source_dir):
        for f in files:
            fp = os.path.join(root, f)
            total_size += os.path.getsize(fp)
            file_count += 1

    print(f"Files to pack: {file_count}")
    print(f"Total size: {total_size / (1024**3):.2f} GB")
    print()
    print("Packing... (this may take a few minutes)")

    # Create tar.gz archive
    with tarfile.open(archive_path, "w:gz") as tar:
        # Add with relative path from project root
        for item in source_dir.iterdir():
            arcname = f"data/output/nplike_v2/{item.name}"
            print(f"  Adding: {arcname}")
            tar.add(item, arcname=arcname)

    # Report result
    archive_size = archive_path.stat().st_size
    print()
    print(f"=" * 60)
    print(f"Done!")
    print(f"Archive: {archive_path}")
    print(f"Size: {archive_size / (1024**3):.2f} GB")
    print(f"Compression ratio: {archive_size / total_size * 100:.1f}%")
    print(f"=" * 60)
    print()
    print("Next steps:")
    print(f"  1. Transfer to server: scp {archive_name} user@server:/path/to/")
    print(f"  2. On server, extract: tar -xzvf {archive_name}")

if __name__ == "__main__":
    main()
