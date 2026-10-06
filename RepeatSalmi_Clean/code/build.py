"""Build only the C++ programs used for Figures 1, 3 and 8."""
from pathlib import Path
import argparse
import shutil
import subprocess

ROOT = Path(__file__).resolve().parent
EXPERIMENT = ROOT / 'experiments/paper_figure3'
BASE = 'numerics atmosphere beam stopping radiation compton redistribution redistribution_bins checkpoint gaunt snapshot anderson'.split()
ATMOSPHERE = 'main kernel_library hydrostatic heating opacity transfer temperature balanced_temperature coupled_probe atlas_numerics'.split()
SOURCES = {
    'atmosphere': [EXPERIMENT / (n + '.cpp') for n in ATMOSPHERE] + [ROOT / 'src' / (n + '.cpp') for n in BASE],
    'figure1': [ROOT / 'app/figure1.cpp'] + [ROOT / 'src' / (n + '.cpp') for n in ('stopping', 'atmosphere', 'numerics')],
}

def build(target, output, compiler=None):
    compiler = compiler or shutil.which('g++')
    if not compiler:
        candidate = Path(r'C:\Program Files (x86)\Dev-Cpp\MinGW64\bin\g++.exe')
        if candidate.is_file():
            compiler = str(candidate)
    if not compiler:
        raise RuntimeError('Install GCC with C++14 and OpenMP, or supply --compiler.')
    output = Path(output).resolve()
    output.parent.mkdir(parents=True, exist_ok=True)
    command = [compiler, '-std=c++14', '-O2', '-fopenmp', '-I' + str(ROOT / 'include'), '-I' + str(EXPERIMENT)]
    if Path(compiler).suffix.lower() == '.exe':
        command += ['-static']
    subprocess.run(command + [str(f) for f in SOURCES[target]] + ['-o', str(output)], check=True)
    print('Built', output)
    return output

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--target', choices=SOURCES, default='atmosphere')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--compiler')
    args = parser.parse_args()
    import os
    suffix = '.exe' if os.name == 'nt' else ''
    build(args.target, args.output or ROOT / 'build' / ('salmi_' + args.target + suffix), args.compiler)
