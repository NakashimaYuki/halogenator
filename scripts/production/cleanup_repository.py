#!/usr/bin/env python3
"""
Repository cleanup and reorganization script.
"""

import shutil
from pathlib import Path
import json

class RepositoryCleanup:
    def __init__(self, repo_root):
        self.repo_root = Path(repo_root)
        self.dry_run = True  # Safety first
        self.actions = []

    def analyze(self):
        """Analyze current state and plan actions."""

        print("="*80)
        print("REPOSITORY CLEANUP ANALYSIS")
        print("="*80)

        # Find all .md files (exclude README.md and CHANGELOG.md)
        exclude_md = {'README.md', 'CHANGELOG.md'}
        md_files = [f for f in self.repo_root.glob("*.md") if f.name not in exclude_md]
        print(f"\n[MD FILES] Markdown files in root: {len(md_files)} (excluding {', '.join(exclude_md)})")
        for f in md_files:
            self.actions.append({
                'type': 'move',
                'source': str(f),
                'dest': f"docs/{f.name}",
                'reason': 'Documentation to docs/'
            })

        # Find test directories/files to delete
        test_patterns = [
            "data/output/transforms/OPT_*",
            "data/output/transforms/TEST_*",
            "data/output/transforms/VALIDATION_*",
            "data/test/",
            "data/viz/",
            "data/viz_v2/",
            "data/viz_base_libs/",
            "tmp/"
        ]

        for pattern in test_patterns:
            matches = list(self.repo_root.glob(pattern))
            for match in matches:
                if match.exists():
                    size = self.get_size(match)
                    self.actions.append({
                        'type': 'delete',
                        'path': str(match),
                        'size': size,
                        'reason': 'Test/temporary data'
                    })

        # Find .log files
        log_files = list(self.repo_root.glob("*.log"))
        print(f"\n[LOG FILES] Log files in root: {len(log_files)}")
        for f in log_files:
            self.actions.append({
                'type': 'move',
                'source': str(f),
                'dest': f"logs/{f.name}",
                'reason': 'Log to logs/'
            })

        # Find standalone .py scripts
        py_files = [
            f for f in self.repo_root.glob("*.py")
            if f.name not in ['setup.py', '__init__.py']
        ]
        print(f"\n[PY FILES] Python scripts in root: {len(py_files)}")
        for f in py_files:
            # Categorize
            if 'test' in f.name.lower() or 'validate' in f.name.lower():
                dest = f"scripts/archive/{f.name}"
            elif 'diagnose' in f.name.lower() or 'analyze' in f.name.lower():
                dest = f"scripts/diagnosis/{f.name}"
            elif 'optimize' in f.name.lower() or 'batch' in f.name.lower():
                dest = f"scripts/production/{f.name}"
            else:
                dest = f"scripts/archive/{f.name}"

            self.actions.append({
                'type': 'move',
                'source': str(f),
                'dest': dest,
                'reason': 'Script organization'
            })

        # Summary
        print(f"\n" + "="*80)
        print(f"PLANNED ACTIONS SUMMARY")
        print(f"="*80)

        action_counts = {}
        for action in self.actions:
            action_type = action['type']
            action_counts[action_type] = action_counts.get(action_type, 0) + 1

        for action_type, count in action_counts.items():
            print(f"  {action_type.upper()}: {count} items")

        # Calculate space to be freed
        delete_size = sum(
            action.get('size', 0)
            for action in self.actions
            if action['type'] == 'delete'
        )
        print(f"\n[DISK SPACE] Space to be freed: {delete_size / (1024**3):.2f} GB")

        return self.actions

    def get_size(self, path):
        """Get size of file or directory in bytes."""
        path = Path(path)
        if path.is_file():
            return path.stat().st_size
        elif path.is_dir():
            return sum(f.stat().st_size for f in path.rglob('*') if f.is_file())
        return 0

    def execute(self, dry_run=True):
        """Execute planned actions."""

        self.dry_run = dry_run

        mode_str = "DRY RUN" if dry_run else "EXECUTING"
        print(f"\n{'='*80}")
        print(f"{mode_str} - Repository Cleanup")
        print(f"{'='*80}")

        success_count = 0
        error_count = 0

        for i, action in enumerate(self.actions, 1):
            print(f"\n[{i}/{len(self.actions)}] {action['type'].upper()}: {action.get('source', action.get('path'))}")

            try:
                if action['type'] == 'move':
                    if not dry_run:
                        source = Path(action['source'])
                        dest = self.repo_root / action['dest']
                        dest.parent.mkdir(parents=True, exist_ok=True)
                        shutil.move(str(source), str(dest))
                    print(f"  -> {action['dest']}")
                    success_count += 1

                elif action['type'] == 'delete':
                    path = Path(action['path'])
                    size_mb = action.get('size', 0) / (1024**2)
                    if not dry_run:
                        if path.is_dir():
                            shutil.rmtree(path)
                        else:
                            path.unlink()
                    print(f"  X Deleted ({size_mb:.1f} MB)")
                    success_count += 1

            except Exception as e:
                print(f"  X ERROR: {e}")
                error_count += 1

        # Summary
        print(f"\n{'='*80}")
        print(f"CLEANUP {'DRY RUN' if dry_run else 'EXECUTION'} COMPLETE")
        print(f"{'='*80}")
        print(f"  Success: {success_count}/{len(self.actions)}")
        print(f"  Errors: {error_count}")

        if dry_run:
            print(f"\n[WARNING] This was a DRY RUN. No files were actually moved or deleted.")
            print(f"   Review the actions above. To execute for real, run:")
            print(f"   python scripts/production/cleanup_repository.py --execute")
        else:
            print(f"\n[OK] Repository cleanup complete!")

            # Save cleanup report
            report_file = self.repo_root / "docs" / "reports" / "cleanup_report.json"
            report_file.parent.mkdir(parents=True, exist_ok=True)
            with open(report_file, 'w') as f:
                json.dump({
                    'actions': self.actions,
                    'success_count': success_count,
                    'error_count': error_count
                }, f, indent=2)
            print(f"   Report saved to: {report_file}")

def main():
    import argparse
    parser = argparse.ArgumentParser(description='Clean up and reorganize repository')
    parser.add_argument('--execute', action='store_true', help='Execute cleanup (default is dry run)')
    parser.add_argument('--repo', default='.', help='Repository root path')
    args = parser.parse_args()

    cleanup = RepositoryCleanup(args.repo)
    cleanup.analyze()
    cleanup.execute(dry_run=not args.execute)

if __name__ == '__main__':
    main()
