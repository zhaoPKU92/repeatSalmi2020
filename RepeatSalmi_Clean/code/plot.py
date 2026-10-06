"""Plot computed Figures 1, 3 and 8 without bundled results or digitized curves."""
from pathlib import Path
import argparse
import json
import numpy as np
from scipy.interpolate import CubicSpline
from numpy.polynomial.legendre import leggauss
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ENERGIES = [.1, .3, 1., 3., 10.]
COLORS = ['black', '#1559e8', '#168b21', '#ed8a13', '#e82424']

def csv(path):
    return np.atleast_1d(np.genfromtxt(path, delimiter=',', names=True))

def atmosphere(directory):
    status = json.loads((directory / 'status.json').read_text())
    if not status.get('discrete_converged'):
        raise ValueError('Atmosphere has not passed equilibrium tolerances: ' + str(directory))
    data = csv(directory / 'intensity_grid.csv')
    energy = np.unique(data['energy_keV'])
    mu = np.unique(data['mu'])
    data = np.sort(data, order=['energy_keV', 'mu'])
    intensity = data['I_E'].reshape(len(energy), len(mu))
    if not np.all(np.isfinite(intensity) & (intensity > 0)):
        raise ValueError('Invalid outgoing intensities')
    nodes, weights = leggauss(len(mu))
    if not np.allclose(mu, (nodes + 1) / 2, rtol=1e-12, atol=1e-14):
        raise ValueError('Unexpected angular quadrature')
    flux = 2 * np.pi * (intensity @ (mu * weights / 2))
    return status, energy, mu, intensity, flux

def figure1(directory, destination):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    for gamma, color in zip([2, 10, 50, 500], COLORS[1:]):
        for brems in [0, 1]:
            data = csv(directory / ('g%d_b%d.csv' % (gamma, brems)))
            style = '--' if brems else '-'
            axes[0].loglog(data['tau'], data['phi'], style, color=color, label=r'$\gamma_0=%d$' % gamma if not brems else None)
            axes[1].loglog(data['tau'], data['gamma'], style, color=color)
    axes[0].set(xlabel=r'$\tau$', ylabel=r'$\Phi=-d\gamma/d\tau$')
    axes[1].set(xlabel=r'$\tau$', ylabel=r'$\gamma$')
    axes[0].legend()
    fig.suptitle('Grey atmosphere: solid Coulomb/Langmuir, dashed with bremsstrahlung')
    fig.tight_layout()
    fig.savefig(destination / 'figure1.pdf')
    plt.close(fig)

def figure3(work, destination):
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.5))
    for name, label in [('nonheated', 'Non-heated'), ('gamma500', r'$\gamma_0=500$'), ('gamma50', r'$\gamma_0=50$'), ('gamma10', r'$\gamma_0=10$')]:
        directory = work / (name + '_a11')
        _, energy, _, _, flux = atmosphere(directory)
        structure = csv(directory / 'structure.csv')
        axes[0].loglog(structure['mass'], structure['T'], label=label)
        axes[1].loglog(energy, energy * flux, label=label)
    axes[0].set(xlabel=r'$m\;(\mathrm{g\,cm^{-2}})$', ylabel=r'$T\;(\mathrm{K})$')
    axes[1].set(xlabel=r'$E\;(\mathrm{keV})$', ylabel=r'$E F_E\;(\mathrm{erg\,cm^{-2}\,s^{-1}})$', xlim=(.01, 100))
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(destination / 'figure3.pdf')
    plt.close(fig)

def beaming(directory):
    status, energy, mu, intensity, _ = atmosphere(directory)
    if len(mu) != 11:
        raise ValueError('Use the final 11-angle atmosphere for beaming.')
    if energy[0] > min(ENERGIES) or energy[-1] < max(ENERGIES):
        raise ValueError('Requested photon energies are outside the model grid.')
    sampled = np.exp(CubicSpline(np.log(energy), np.log(intensity), axis=0)(np.log(ENERGIES)))
    # The positive Gauss nodes stop below 1. Extend the angular spline to
    # normal emission; do not invent intensities below the smallest mu node.
    normal = CubicSpline(mu, sampled, axis=1)(1.)
    if np.any(normal <= 0):
        raise ValueError('Invalid normal-emission extrapolation')
    return status, mu, sampled / normal[:, None]

def figure8(work, destination, deltas):
    cold_status, cold_mu, cold = beaming(work / 'nonheated_a11')
    for delta in deltas:
        status, mu, heated = beaming(work / ('delta%d_a11' % delta))
        if status['delta'] != delta or cold_status['delta'] != 0:
            raise ValueError('Unexpected model identity')
        fig = plt.figure(figsize=(10, 4.5))
        polar = fig.add_subplot(121, projection='polar')
        cart = fig.add_subplot(122)
        polar.set_theta_zero_location('N')
        polar.set_thetamin(0)
        polar.set_thetamax(90)
        for j, (energy, color) in enumerate(zip(ENERGIES, COLORS)):
            for x, values, style in [(cold_mu, cold[j], '--'), (mu, heated[j], '-')]:
                knots = np.r_[x, 1.]
                ratio = np.r_[values, 1.]
                fine = np.linspace(x[0], 1., 500)
                curve = CubicSpline(knots, ratio)(fine)
                if np.any(curve <= 0):
                    raise ValueError('Angular interpolation became nonpositive')
                cart.semilogy(fine, curve, style, color=color, label='%g keV' % energy if style == '-' else None)
                polar.plot(np.arccos(fine), np.log10(curve), style, color=color)
        cart.set(xlabel=r'$\mu=\cos\theta$', ylabel=r'$I_E(\mu)/I_E(1)$', xlim=(0, 1))
        cart.legend()
        polar.set_title(r'$\log_{10}[I_E(\mu)/I_E(1)]$')
        fig.suptitle(r'$\delta=%d$: solid heated, dashed non-heated' % delta)
        fig.tight_layout()
        fig.savefig(destination / ('%d.pdf' % delta))
        plt.close(fig)

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('figure', type=int, choices=[1, 3, 8])
    parser.add_argument('--work', type=Path, default=Path(__file__).resolve().parent / 'work')
    parser.add_argument('--output', type=Path)
    parser.add_argument('--delta', type=int, choices=[1, 2, 3])
    args = parser.parse_args()
    destination = args.output or args.work / 'figures'
    destination.mkdir(parents=True, exist_ok=True)
    if args.figure == 1:
        figure1(args.work / 'figure1', destination)
    elif args.figure == 3:
        figure3(args.work, destination)
    else:
        figure8(args.work, destination, [args.delta] if args.delta else [1, 2, 3])
    print('Saved plots in', destination)

if __name__ == '__main__':
    main()
