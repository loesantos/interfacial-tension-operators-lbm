#!/usr/bin/env python3
"""LaTeX tables for the paper from the Lotes/*.csv files and the translating-drop derivation."""
import csv, json, math, sys
from pathlib import Path
L = Path(sys.argv[1]); DER = json.load(open(sys.argv[2]))
OUT = Path(__file__).resolve().parent
OPS = ['P1', 'P3 chi=-0.4', 'P3 chi=0.5', 'P3 chi=2', 'F1', 'F2', 'F3']
LAB = {'P1': 'PG', 'P3 chi=-0.4': 'PL, RP', 'P3 chi=0.5': 'PL, LKR', 'P3 chi=2': 'PL, LVK',
       'F1': 'F1', 'F2': 'F2', 'F3': 'F3'}
rd = lambda f: list(csv.DictReader(open(L / f)))
num = float

def e(x, d=1):
    if x == 0: return '0'
    s = f'{x:.{d}e}'; m, ex = s.split('e'); ex = int(ex)
    return f'${m}\\times10^{{{ex}}}$'

# ---- E1 -------------------------------------------------------------------------------------------------------
def fit(Rs, sg):
    import numpy as np
    X = np.array([[1, 1 / R, 1 / R ** 2] for R in Rs]); y = np.array(sg)
    a, *_ = np.linalg.lstsq(X, y, rcond=None)
    pr = X @ a; r2 = 1 - ((y - pr) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return a[0], a[1] / a[0], a[2] / a[0], r2
E1 = rd('E1.csv')
t = ['\\begin{tabular}{lrrrrrrr}', '\\toprule',
     'Operator & $R=10$ & $R=25$ & $R=50$ & $\\sigma_\\infty$ & $\\sigma_\\infty/\\sigma_0-1$ & $\\delta$ & $c_2$\\\\', '\\midrule']
for o in OPS:
    cs = sorted([c for c in E1 if c['rotulo'] == o], key=lambda c: num(c['R_ini']))
    by = {int(num(c['R_ini'])): num(c['sigma_meiafinal']) for c in cs}
    s0, d, c2, r2 = fit([num(c['R_meiafinal']) for c in cs], [num(c['sigma_meiafinal']) for c in cs])
    t.append(f"{LAB[o]} & {by[10]:.5f} & {by[25]:.5f} & {by[50]:.5f} & {s0:.5f} & ${100*(s0/0.1-1):+.2f}\\,\\%$ & ${d:+.3f}$ & ${c2:.2f}$\\\\")
t += ['\\bottomrule', '\\end{tabular}']
(OUT / 'tab_E1.tex').write_text('\n'.join(t) + '\n')

# ---- E4 -------------------------------------------------------------------------------------------------------
E4 = rd('E4.csv')
taus = sorted(set(num(r['tau']) for r in E4))
t = ['\\begin{tabular}{l' + 'c' * len(OPS) + '}', '\\toprule',
     '$\\tau$ & ' + ' & '.join(LAB[o] for o in OPS) + '\\\\', '\\midrule']
for tau in taus:
    cel = []
    for o in OPS:
        rs = [r for r in E4 if r['rotulo'] == o and num(r['tau']) == tau]
        ok = [num(r['sigma_nominal']) for r in rs if r['estado'] == 'estavel' and num(r['neg_total_max']) == 0]
        st = [num(r['sigma_nominal']) for r in rs if r['estado'] == 'estavel']
        cel.append(f"{max(ok) if ok else 0:g}\\,/\\,{max(st) if st else 0:g}")
    t.append(f'{tau:g} & ' + ' & '.join(cel) + '\\\\')
t += ['\\bottomrule', '\\end{tabular}']
(OUT / 'tab_E4.tex').write_text('\n'.join(t) + '\n')

# ---- E4v ------------------------------------------------------------------------------------------------------
E2 = rd('E2.csv')
stat = {r['rotulo']: num(r['umax_sigma']) for r in E2 if num(r['tau']) == 1 and num(r['sigma_nominal']) == 0.1}
# slip from galilean.dat files
slip, umed = {}, {}
for o in OPS:
    for v in ('0.1', '0.25'):
        rows = [list(map(float, l.split())) for l in open(L / 'E4v' / o / f'v0={v}' / 'galilean.dat') if not l.startswith('#')]
        h = rows[len(rows) // 2:]
        slip[(o, v)] = sum(r[13] - r[10] for r in h) / len(h) / float(v)
        umed[(o, v)] = sum(r[10] for r in h) / len(h) / float(v)       # domain-mean velocity / v0
t = ['\\begin{tabular}{lrrrrrrrr}', '\\toprule',
     '& \\multicolumn{3}{c}{$|\\uvec-\\langle\\uvec\\rangle|_{\\max}/\\sigma$} & $\\sum\\bm F_x$ & \\multicolumn{2}{c}{$\\langle u_x\\rangle/v_0$} & \\multicolumn{2}{c}{slip$/v_0$}\\\\',
     '\\cmidrule(lr){2-4}\\cmidrule(lr){5-5}\\cmidrule(lr){6-7}\\cmidrule(lr){8-9}',
     'Operator & $v_0=0$ & $0.1$ & $0.25$ & $v_0=0.1$ & $v_0=0.1$ & $0.25$ & $v_0=0.1$ & $0.25$\\\\', '\\midrule']
for o in OPS:
    d = {x[1]: x for x in DER if x[0] == o}
    sf = d['0.1'][3]
    sfs = '$\\approx0$' if abs(sf) < 1e-12 else e(sf)
    t.append(f"{LAB[o]} & {e(stat[o])} & {e(d['0.1'][2])} & {e(d['0.25'][2])} & {sfs} & {umed[(o,'0.1')]:.3f} & {umed[(o,'0.25')]:.3f} & {e(slip[(o,'0.1')])} & {e(slip[(o,'0.25')])}\\\\")
t += ['\\bottomrule', '\\end{tabular}']
(OUT / 'tab_E4v.tex').write_text('\n'.join(t) + '\n')

# ---- E3 widths ------------------------------------------------------------------------------------------------
E3 = rd('E3.csv')
th = sorted(set(num(r['theta']) for r in E3))
t = ['\\begin{tabular}{l' + 'r' * len(th) + 'r}', '\\toprule',
     'Operator & ' + ' & '.join(f'${x:.1f}^\\circ$' for x in th) + ' & $\\max|\\sigma_{\\mathrm{pop}}/\\sigma-1|$\\\\', '\\midrule']
for o in OPS:
    rs = sorted([r for r in E3 if r['rotulo'] == o], key=lambda r: num(r['theta']))
    w = [2 * math.log(9) / num(r['K_med']) for r in rs]
    sp = max(abs(num(r['sigma_pop_med']) / num(r['sigma_nominal']) - 1) for r in rs) if o.startswith('P') else None
    t.append(f"{LAB[o]} & " + ' & '.join(f'{x:.3f}' for x in w) + ' & ' + (e(sp) if sp else '---') + '\\\\')
t += ['\\bottomrule', '\\end{tabular}']
(OUT / 'tab_E3.tex').write_text('\n'.join(t) + '\n')
print('ok')
