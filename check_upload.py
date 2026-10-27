#!/usr/bin/env python3
"""Check upload contents without executing or modifying the scientific analyses."""
from pathlib import Path
import argparse
import ast
import csv
import hashlib
import importlib.metadata
import json
import re
import shutil
import subprocess
import sys
import tempfile

ROOT = Path(__file__).resolve().parent

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--environment', action='store_true', help='Also check installed R/Python prerequisites; does not install anything.')
    parser.add_argument('--profile', choices=['standard', 'extended'], default='standard')
    parser.add_argument('--rscript', default=shutil.which('Rscript'))
    parser.add_argument('--report', type=Path, help='Optional JSON report; use a location outside the release folder.')
    args = parser.parse_args()
    errors = []
    required = ['README.md', 'docs/RELEASE_CHECKLIST.md', 'docs/ARCHIVED_RECORDS.md', 'LICENSE_STATUS.md',
                'CITATION.cff', '.gitignore', 'reproduce.py', 'verify_package.py', 'renv.lock',
                'environment/requirements.txt', 'environment/requirements-extended.txt',
                'environment/restore.R', 'metadata/FIGURE_TABLE_CODE_MAP.csv',
                'metadata/CODE_COMPLETENESS.csv', 'docs/CODE_COMPLETENESS.md',
                'docs/DATA_PROVENANCE.md', 'docs/VALIDATION.md',
                'licenses/THIRD_PARTY_NOTICES.md', 'manifest.csv', 'SHA256SUMS']
    for relative in required:
        if not (ROOT / relative).is_file():
            errors.append({'path': relative, 'problem': 'missing required release file'})

    entries = []
    tree = ast.parse((ROOT / 'reproduce.py').read_text(encoding='utf-8-sig'))
    for node in ast.walk(tree):
        if isinstance(node, ast.Call) and isinstance(node.func, ast.Name) and node.func.id == 'stage':
            if len(node.args) >= 2 and all(isinstance(value, ast.Constant) for value in node.args[:2]):
                name, script = [value.value for value in node.args[:2]]
                relative = script if script.startswith('code/') else 'code/project/' + script
                entries.append({'step': name, 'script': relative, 'present': (ROOT / relative).is_file()})
    for name in ['validate_results.py', 'export_results.py']:
        relative = 'code/' + name
        entries.append({'step': name.removesuffix('.py'), 'script': relative, 'present': (ROOT / relative).is_file()})
    for entry in entries:
        if not entry['present']:
            errors.append({'path': entry['script'], 'problem': 'missing execution entry'})

    file_count = 0
    largest = {'path': None, 'bytes': 0}
    manifest_path = ROOT / 'manifest.csv'
    if manifest_path.is_file():
        with manifest_path.open(encoding='utf-8-sig', newline='') as stream:
            for row in csv.DictReader(stream):
                file_count += 1
                path = (ROOT / row['path']).resolve()
                if ROOT not in path.parents or not path.is_file():
                    errors.append({'path': row['path'], 'problem': 'missing or unsafe manifest entry'})
                    continue
                if path.stat().st_size != int(row['bytes']) or digest(path) != row['sha256']:
                    errors.append({'path': row['path'], 'problem': 'size or SHA256 mismatch'})
                if path.stat().st_size > largest['bytes']:
                    largest = {'path': row['path'], 'bytes': path.stat().st_size}
                if path.stat().st_size > 100 * 1024 ** 2:
                    errors.append({'path': row['path'], 'problem': 'exceeds 100 MiB regular Git file limit'})
        sums = ROOT / 'SHA256SUMS'
        if sums.is_file():
            expected = [line.split('  ', 1)[0] for line in sums.read_text(encoding='utf-8').splitlines()
                        if line.endswith('  manifest.csv')]
            if len(expected) != 1 or digest(manifest_path) != expected[0]:
                errors.append({'path': 'manifest.csv', 'problem': 'manifest checksum mismatch'})

    syntax_count = 0
    for path in [*ROOT.glob('*.py'), *(ROOT / 'code').rglob('*.py')]:
        try:
            ast.parse(path.read_text(encoding='utf-8-sig'))
            syntax_count += 1
        except SyntaxError as exc:
            errors.append({'path': path.relative_to(ROOT).as_posix(), 'problem': str(exc.msg)})

    environment = {'checked': args.environment, 'profile': args.profile}
    if args.environment:
        needed = {'numpy': '2.3.5', 'pandas': '3.0.1'}
        if args.profile == 'extended':
            needed['gseapy'] = '1.3.1'
        environment['python'] = {}
        for package, expected in needed.items():
            try:
                actual = importlib.metadata.version(package)
            except importlib.metadata.PackageNotFoundError:
                actual = None
            environment['python'][package] = {'required': expected, 'installed': actual}
            if actual != expected:
                errors.append({'path': package, 'problem': 'Python dependency version mismatch'})
        if not args.rscript:
            errors.append({'path': 'Rscript', 'problem': 'not found; supply --rscript'})
        else:
            import os
            env = os.environ.copy()
            if (ROOT / 'environment/r_library').is_dir():
                env['R_LIBS_USER'] = str(ROOT / 'environment/r_library')
            with tempfile.TemporaryDirectory(prefix='sci_prerequisites_') as directory:
                check = subprocess.run([args.rscript, '--vanilla', str(ROOT / 'code/check_environment.R'), directory],
                                       capture_output=True, text=True, encoding='utf-8', errors='replace', env=env)
                environment['R_prerequisites_passed'] = check.returncode == 0
                if check.returncode:
                    errors.append({'path': 'R environment', 'problem': (check.stderr or check.stdout)[-2000:]})

    report = {'package_check_passed': not errors, 'execution_entries_checked': len(entries),
              'missing_execution_entries': sum(not entry['present'] for entry in entries),
              'manifest_files_checked': file_count, 'python_sources_parsed': syntax_count,
              'largest_release_file': largest, 'execution_entries': entries,
              'environment': environment, 'errors': errors,
              'original_code_license': 'PENDING_AUTHOR_DECISION', 'publication_status': 'PUBLIC_GITHUB_ACCESS_VERIFIED_2026_10_27',
              'availability_note': 'Recorded public-access verification; this local integrity check does not query the remote repository.',
              'scope': 'Package integrity and declared entry points; not a new complete numerical rerun or clean-machine installation test.'}
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding='utf-8')
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report['package_check_passed'] else 1

if __name__ == '__main__':
    raise SystemExit(main())
