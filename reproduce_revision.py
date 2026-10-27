#!/usr/bin/env python3
"""Run the October revision extensions in a new, separate output directory."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess
import sys

ROOT = Path(__file__).resolve().parent


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--task', choices=['verify', 'pilot', 'full', 'figures'], default='verify')
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--rscript', default=shutil.which('Rscript'))
    args = p.parse_args()
    out = args.output.resolve()
    if out.exists():
        p.error('Use a new output directory; existing files will not be overwritten.')
    if out == ROOT or out in ROOT.parents or (ROOT in out.parents and ROOT/'outputs' not in out.parents):
        p.error('Inside the repository, outputs must be below outputs/.')
    if args.task != 'verify' and not args.rscript:
        p.error('Rscript is required; supply --rscript.')
    out.mkdir(parents=True)
    work = out/'work'
    logs = out/'logs'
    logs.mkdir()
    env = os.environ.copy()
    env['PYTHONIOENCODING'] = 'utf-8'
    if (ROOT/'environment/r_library').is_dir():
        env['R_LIBS_USER'] = str(ROOT/'environment/r_library')
    commands = []

    def run(name, command):
        print(f'Running {name}', flush=True)
        with (logs/f'{name}.log').open('w', encoding='utf-8') as stream:
            result = subprocess.run(command, cwd=work, env=env, stdout=stream,
                                    stderr=subprocess.STDOUT, text=True)
        commands.append({'step': name, 'returncode': result.returncode})
        if result.returncode:
            raise RuntimeError(f'{name} failed; see {logs/name}.log')

    def r(name, *extra):
        run(Path(name).stem, [args.rscript, '--vanilla', str(work/'scripts'/name), str(work), *extra])

    if args.task == 'figures':
        src = ROOT/'extensions/figure_redraw'
        work.mkdir()
        shutil.copytree(src/'source_data', work/'source_data')
        shutil.copy2(src/'make_figures.R', work/'make_figures.R')
        run('make_figures', [args.rscript, '--vanilla', str(work/'make_figures.R'), str(work)])
        detail = {'figures': len(list(work.glob('Figure*.pdf'))),
                  'tiff_files': len(list(work.glob('Figure*.tif'))),
                  'scientific_models_refitted': False}
        if detail['figures'] != 7 or detail['tiff_files'] != 7:
            raise RuntimeError('Expected seven PDF and seven TIFF outputs.')
    else:
        src = ROOT/'extensions/nested_validation'
        if args.task in {'verify', 'pilot'}:
            shutil.copytree(src, work)
        else:
            work.mkdir()
            for folder in ['inputs', 'scripts']:
                shutil.copytree(src/folder, work/folder)
            for name in ['PROTOCOL.md', 'requirements-analysis.txt']:
                shutil.copy2(src/name, work/name)
            for folder in ['audit', 'tables', 'logs', 'figures']:
                (work/folder).mkdir(exist_ok=True)
            shutil.copy2(src/'audit/protocol_input_lock.json', work/'audit/protocol_input_lock.json')
        if args.task in {'pilot', 'full'}:
            r('run_nested.R', 'pilot')
        if args.task == 'full':
            r('run_nested.R', 'observed')
            r('run_nested.R', 'permutations')
            r('run_clinical_resampling.R')
            r('summarize_results.R')
            r('make_figures.R')
        run('verify_results', [sys.executable, str(work/'scripts/verify_results.py')])
        detail = json.loads((work/'audit/analysis_verification.json').read_text(encoding='utf-8'))
    report = {'passed': True, 'task': args.task, 'commands': commands,
              'details': detail, 'source_reference_files_modified': False,
              'full_nested_and_resampling_rerun': args.task == 'full'}
    (out/'validation.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(json.dumps({'passed': True, 'task': args.task, 'report': str(out/'validation.json')}))


if __name__ == '__main__':
    main()
