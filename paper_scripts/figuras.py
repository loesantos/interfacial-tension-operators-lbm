#!/usr/bin/env python3
"""Figures for the interfacial-tension-operator paper, read from the Lotes/*.csv files.

    python3 figuras.py <path to Lotes>        # writes figures/*.pdf next to this script
"""
import csv, json, math, sys
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

LOTES = Path(sys.argv[1]) if len(sys.argv) > 1 else Path('Lotes')
OUT = Path(__file__).resolve().parent / 'figures'
OUT.mkdir(exist_ok=True)

plt.rcParams.update({
    'font.family': 'serif', 'font.size': 9, 'axes.labelsize': 9, 'legend.fontsize': 7.5,
    'xtick.labelsize': 8, 'ytick.labelsize': 8, 'axes.linewidth': 0.6, 'lines.linewidth': 1.4,
    'xtick.major.width': 0.6, 'ytick.major.width': 0.6, 'axes.edgecolor': '#555555',
    'axes.grid': True, 'grid.color': '#e4e4e0', 'grid.linewidth': 0.5, 'savefig.bbox': 'tight',
    'mathtext.fontset': 'cm'})

#  Categorical slots in fixed order ( validated: adjacent CVD dE >= 9.1 ), plus marker and dash as
#  secondary encoding.  Perturbations dashed, forces solid.
OPS = [
    ('P1',          r'PG ($\lambda=0$)',            '#2a78d6', 'o', '--'),
    ('P3 chi=-0.4', r'PL, $\chi=-2/5$ (RP)',        '#eb6834', '^', '--'),
    ('P3 chi=0.5',  r'PL, $\chi=1/2$ (LKR)',        '#1baf7a', 'v', '--'),
    ('P3 chi=2',    r'PL, $\chi=2$ (LVK)',          '#eda100', 'D', '--'),
    ('F1',          'F1 (CSF)',                     '#e87ba4', 'o', '-'),
    ('F2',          r'F2 ($\nabla\cdot\mathsf{T}$)', '#008300', 's', '-'),
    ('F3',          r'F3 ($\mu\nabla\phi$)',        '#4a3aa7', '^', '-'),
]
INK = '#2b2b2b'


def le(nome):
    with open(LOTES / nome) as f:
        return [{k: (float(v) if k not in ('pasta', 'rotulo', 'operador', 'estado') else v)
                 for k, v in r.items() if v not in (None, '')} for r in csv.DictReader(f)]


def estilo(ax):
    ax.tick_params(colors=INK)
    for s in ('top', 'right'):
        ax.spines[s].set_visible(False)


#------------------------------------------------------------------------------------------------------------------#
# Fig. 1: spurious velocity against chi ( drop, tau = 1 and 0.6; inclined flat interface )
#------------------------------------------------------------------------------------------------------------------#

e2chi, e3chi = le('E2chi.csv'), le('E3chi.csv')
fig, ax = plt.subplots(figsize=(3.4, 2.6))
series = [
    ([r for r in e2chi if r['tau'] == 1.0], 'drop, $\\tau=1$', '#2a78d6', 'o', '-'),
    ([r for r in e2chi if r['tau'] == 0.6], 'drop, $\\tau=0.6$', '#eb6834', 's', '-'),
    (e3chi, 'flat, $\\theta=26.6^\\circ$, $\\tau=1$', '#1baf7a', '^', '--'),
]
for rows, lab, c, m, ls in series:
    rows = sorted(rows, key=lambda r: r['chi'])
    y = [r['umax_sigma'] if 'umax_sigma' in r else r['umax/sigma_nominal'] for r in rows]
    ax.plot([r['chi'] for r in rows], y, ls, color=c, marker=m, ms=4, label=lab, mec='white', mew=0.6)
for chi, nome in ((-0.4, 'RP'), (0.5, 'LKR'), (2.0, 'LVK')):
    ax.axvline(chi, color='#9a9a94', lw=0.6, ls=':', zorder=0)
    ax.text(chi, 0.2, nome, ha='center', va='bottom', fontsize=7, color='#666660',
            bbox=dict(fc='white', ec='none', pad=0.5))
ax.set_yscale('log')
ax.set_xlabel(r'$\chi$')
ax.set_ylabel(r'$|\mathbf{u}|_{\max}/\sigma$')
ax.set_ylim(8e-4, 0.3)
ax.legend(frameon=False, loc='lower right')
estilo(ax)
fig.savefig(OUT / 'fig_chi.pdf')
plt.close(fig)

#------------------------------------------------------------------------------------------------------------------#
# Fig. 2: spurious velocity at a straight inclined interface against theta
#------------------------------------------------------------------------------------------------------------------#

e3 = le('E3.csv')
fig, ax = plt.subplots(figsize=(3.4, 2.7))
for rot, lab, c, m, ls in OPS:
    rows = sorted([r for r in e3 if r['rotulo'] == rot], key=lambda r: r['theta'])
    y = [max(r['umax/sigma_nominal'], 1e-14) for r in rows]
    ax.plot([r['theta'] for r in rows], y, ls, color=c, marker=m, ms=4, label=lab, mec='white', mew=0.6)
ax.set_yscale('log')
ax.set_xlabel(r'interface normal angle $\theta$ (deg)')
ax.set_ylabel(r'$|\mathbf{u}|_{\max}/\sigma$')
ax.set_xticks([0, 14.04, 18.43, 26.57, 33.69, 45])
ax.set_xticklabels(['0', '14', '18', '27', '34', '45'])
ax.set_ylim(1e-15, 0.2)
ax.legend(frameon=False, ncol=2, loc='upper center', bbox_to_anchor=(0.5, -0.2), columnspacing=1.0)
estilo(ax)
fig.savefig(OUT / 'fig_theta.pdf')
plt.close(fig)

#------------------------------------------------------------------------------------------------------------------#
# Fig. 3: capillary number of the spurious flow, Ca_s = nu |u|max / sigma ( static drop, R = 25 )
#------------------------------------------------------------------------------------------------------------------#

e2 = le('E2.csv')
fig, (a1, a2) = plt.subplots(1, 2, figsize=(6.8, 2.7))
for rot, lab, c, m, ls in OPS:
    rows = sorted([r for r in e2 if r['rotulo'] == rot and r['tau'] == 1.0], key=lambda r: r['sigma_nominal'])
    a1.plot([r['sigma_nominal'] for r in rows], [r['Ca_s'] for r in rows], ls, color=c, marker=m, ms=4,
            label=lab, mec='white', mew=0.6)
    rows = sorted([r for r in e2 if r['rotulo'] == rot and r['sigma_nominal'] == 0.1], key=lambda r: r['tau'])
    a2.plot([r['tau'] for r in rows], [r['Ca_s'] for r in rows], ls, color=c, marker=m, ms=4, mec='white', mew=0.6)
a1.set_xscale('log'); a1.set_yscale('log')
a1.set_xlabel(r'input $\sigma$'); a1.set_ylabel(r'$\mathrm{Ca}_s=\nu|\mathbf{u}|_{\max}/\sigma_{\mathrm{meas}}$')
a1.set_title(r'(a) $\tau=1$', fontsize=9, loc='left')
a2.set_yscale('log')
a2.set_xlabel(r'$\tau$'); a2.set_ylabel(r'$\mathrm{Ca}_s$')
a2.set_title(r'(b) $\sigma=0.1$', fontsize=9, loc='left')
a2.axvline(0.5 + math.sqrt(1 / 12), color='#9a9a94', lw=0.6, ls=':', zorder=0)
a2.text(0.5 + math.sqrt(1 / 12), 3e-2, r'$\Lambda=\frac{1}{12}$', ha='center', fontsize=7, color='#666660',
        bbox=dict(fc='white', ec='none', pad=0.5))
for a in (a1, a2):
    a.set_ylim(5e-6, 6e-2); estilo(a)
fig.legend(*a1.get_legend_handles_labels(), frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.12))
fig.subplots_adjust(wspace=0.3, bottom=0.25)
fig.savefig(OUT / 'fig_ca.pdf')
plt.close(fig)

#------------------------------------------------------------------------------------------------------------------#
# Fig. 4: translating drop -- spurious velocity relative to the mean flow
#------------------------------------------------------------------------------------------------------------------#

deriv = json.load(open(sys.argv[2])) if len(sys.argv) > 2 else []


def historia(rot, v):
    """Time history of a translating-drop run: step, <u_x>/v0, |u-<u>|max/sigma_measured."""
    rows = [list(map(float, l.split())) for l in open(LOTES / 'E4v' / rot / f'v0={v}' / 'galilean.dat')
            if not l.startswith('#')]
    rows = [r for r in rows if r[0] > 0 and abs(r[3]) > 1e-12]
    return [r[0] for r in rows], [r[10] / float(v) for r in rows], [r[7] / r[3] for r in rows]


if deriv:
    stat = {r['rotulo']: r['umax_sigma'] for r in e2 if r['tau'] == 1.0 and r['sigma_nominal'] == 0.1}
    fig, (a1, a2, a3) = plt.subplots(1, 3, figsize=(7.0, 2.6))
    for rot, lab, c, m, ls in OPS:
        pts = sorted([(float(d[1]), d[2]) for d in deriv if d[0] == rot])
        xs = [0.0] + [p[0] for p in pts]
        ys = [stat[rot]] + [p[1] for p in pts]
        a1.plot(xs, ys, ls, color=c, marker=m, ms=4, label=lab, mec='white', mew=0.6)
        t, um, ud = historia(rot, '0.25')
        a3.plot([x / 1e3 for x in t], ud, ls, color=c, lw=1.2)
    #  (b) mean flow: only F1 and F3 depart from v0; the other five operators stay at 1.000
    for rot, c, lsv in (('F3', '#4a3aa7', {'0.1': '-', '0.25': ':'}), ('F1', '#e87ba4', {'0.1': '-', '0.25': ':'})):
        for v in ('0.1', '0.25'):
            t, um, ud = historia(rot, v)
            a2.plot([x / 1e3 for x in t], um, lsv[v], color=c, lw=1.4)
            if not (rot == 'F3' and v == '0.25') and not (rot == 'F1' and v == '0.1'):
                lab = 'F3, both' if rot == 'F3' else f'F1, {v}'
                a2.text(t[-1] / 1e3 + 0.6, um[-1] + (0.03 if v == '0.1' and rot == 'F1' else 0), lab,
                        fontsize=6.5, color=INK, va='bottom' if rot == 'F1' and v == '0.1' else 'center')
    a2.axhline(1.0, color='#9a9a94', lw=0.8, zorder=0)
    a2.text(2, 1.025, r'PG, PL, F2; F1 at $v_0=0.1$', fontsize=6.5, color='#666660', va='bottom')
    a1.set_yscale('log'); a3.set_yscale('log')
    a1.set_xlabel(r'$v_0$'); a1.set_ylabel(r'$|\mathbf{u}-\langle\mathbf{u}\rangle|_{\max}/\sigma$')
    a1.set_title('(a) second-half mean', fontsize=9, loc='left')
    a2.set_xlabel(r'step ($10^3$)'); a2.set_ylabel(r'$\langle u_x\rangle/v_0$')
    a2.set_xlim(0, 40); a2.set_xticks([0, 10, 20, 30]); a2.set_ylim(0, 1.15)
    a2.set_title(r'(b) mean flow, $v_0=0.1$ and $0.25$', fontsize=9, loc='left')
    a3.set_xlabel(r'step ($10^3$)'); a3.set_ylabel(r'$|\mathbf{u}-\langle\mathbf{u}\rangle|_{\max}/\sigma$')
    a3.set_title(r'(c) history, $v_0=0.25$', fontsize=9, loc='left')
    for a in (a1, a2, a3):
        estilo(a)
    fig.legend(*a1.get_legend_handles_labels(), frameon=False, ncol=4, loc='lower center', bbox_to_anchor=(0.5, -0.14))
    fig.subplots_adjust(wspace=0.55, bottom=0.27)
    fig.savefig(OUT / 'fig_translation.pdf')
    plt.close(fig)

print('figures in', OUT)
