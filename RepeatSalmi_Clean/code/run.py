"""Run fresh paper models, or continue a newly computed atmosphere."""
from pathlib import Path
import argparse
import json
import os
import subprocess
from build import ROOT, build
from model_config import MODELS, MASS_GRIDS

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('model', choices=['figure1'] + list(MODELS))
    parser.add_argument('--work', type=Path, default=ROOT / 'work', help='Location for newly generated files')
    parser.add_argument('--angles', type=int, choices=[3, 11], default=3)
    parser.add_argument('--iterations', type=int, help='Default: 160 at 3 angles; 80 at 11 angles')
    parser.add_argument('--threads', type=int, default=2)
    parser.add_argument('--initial', type=Path, help='Newly computed field/checkpoint directory; no bundled results required')
    parser.add_argument('--output', type=Path, help='A new output directory; existing outputs are never overwritten')
    parser.add_argument('--compton', choices=['midpoint', 'moment_v1'], default='midpoint')
    parser.add_argument('--lapack-dll', type=Path, help='Optional Windows LP64 OpenBLAS DLL for the same least-squares solve')
    parser.add_argument('--dry-run', action='store_true', help='Print inputs and command without writing or running')
    args = parser.parse_args()
    if args.threads < 1 or (args.iterations is not None and args.iterations < 1):
        parser.error('Threads and iterations must be positive.')
    work = args.work.resolve()
    target = 'figure1' if args.model == 'figure1' else 'atmosphere'
    binary = work / 'build' / ('salmi_' + target + ('.exe' if os.name == 'nt' else ''))
    output = (args.output or work / (args.model + ('' if target == 'figure1' else '_a' + str(args.angles)))).resolve()
    if output.exists():
        parser.error('Output exists. Choose --output for a new directory, and --initial to continue an earlier calculation.')
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(args.threads), OMP_WAIT_POLICY='PASSIVE', SALMI_PAIR_CACHE_MB='768', SALMI_COMPTON_BIN_INTEGRATION=args.compton)
    # Explicit backend selection prevents an inherited local DLL path from breaking portability.
    env.pop('SALMI_LAPACK_DLL', None)
    if args.lapack_dll:
        dll = args.lapack_dll.resolve()
        if not dll.is_file():
            parser.error('LAPACK DLL does not exist.')
        env['SALMI_LAPACK_DLL'] = str(dll)
    if target == 'figure1':
        command = [str(binary), str(output), '.002', '1e-6']
    else:
        model = MODELS[args.model]
        options = dict(model['options'])
        grid_path = output / 'mass.csv'
        masses = MASS_GRIDS[model['grid']]
        options.update({'--output': str(output).replace('\\', '/'), '--cache': str(work / ('cache_' + args.compton + '_a' + str(args.angles))).replace('\\', '/'),
                        '--mass-grid': str(grid_path).replace('\\', '/'), '--angles': str(args.angles), '--threads': str(args.threads),
                        '--iterations': str(args.iterations or (160 if args.angles == 3 else 80)),
                        '--coefficient-cache-mb': '1024' if args.angles == 3 else '4096',
                        '--temperature-step': '.03' if args.angles == 3 else '.01'})
        if args.initial:
            initial = args.initial.resolve()
            for name in ('structure.csv', 'occupation.bin', 'status.json'):
                if not (initial / name).is_file():
                    parser.error('Initial directory is missing ' + name)
            options.update({'--initial-profile': str(initial / 'structure.csv').replace('\\', '/'), '--initial-field': str(initial).replace('\\', '/')})
        command = [str(binary)] + [item for key, value in options.items() for item in (key, str(value))]
    if args.dry_run:
        print(json.dumps({'command': command, 'initialization': str(args.initial) if args.initial else 'fresh grey/Planck', 'compton': args.compton}, indent=2))
        return
    build(target, binary)
    output.mkdir(parents=True)
    if target == 'atmosphere':
        (work / ('cache_' + args.compton + '_a' + str(args.angles))).mkdir(parents=True, exist_ok=True)
        grid_path.write_text('mass\n' + ''.join(format(m, '.17g') + '\n' for m in masses), encoding='ascii')
    print('Running:', output, flush=True)
    with (output / 'run.log').open('w', encoding='utf-8') as log:
        result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT, env=env)
    print('Exit code:', result.returncode, '; log:', output / 'run.log')
    if target == 'atmosphere' and (output / 'status.json').exists():
        status = json.loads((output / 'status.json').read_text())
        print('Discrete equilibrium:', status['discrete_converged'], '; paper agreement requires separate comparison.')
    raise SystemExit(result.returncode)

if __name__ == '__main__':
    main()
