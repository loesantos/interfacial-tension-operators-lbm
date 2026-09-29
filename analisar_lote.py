#!/usr/bin/env python3
"""Junta os resultados dos lotes de rodar_lote.py em tabelas.

    python3 analisar_lote.py                 # le Lotes/, escreve Lotes/RESULTADOS.md e Lotes/<lote>.csv
    python3 analisar_lote.py --saida Lotes_teste

Para cada caso le resumo.dat ( medias da metade final ) e laplace.dat, de onde tira tambem a
TENDENCIA de |u|max: ( media dos ultimos 20 % ) / ( media entre 50 e 70 % ) - 1.  Acima de 5 % o caso
e marcado como nao convergido ( '*' nas tabelas ).

Tabelas:
    E1      sigma( R ) = sigma_inf ( 1 + delta / R + c2 / R^2 ), por operador
    E2      Ca_s = nu |u|max / sigma, por operador, tau x sigma
    E2chi   |u|max / sigma e a4 contra chi
    E4      mapa de estabilidade tau x sigma:  ok | n ( populacoes totais negativas ) | X ( instavel )
    E4v     gota em translacao: |u - <u>|max / sigma, soma F_x, populacoes de cor negativas
    E3      interface plana inclinada: sigma_pop / sigma, espessura e |u|max / sigma contra theta
    E3chi   idem em theta = 26.6 graus, contra chi
"""

import argparse
import csv
import math
from pathlib import Path

AQUI = Path(__file__).resolve().parent

COLS = ['operador', 'chi', 'xi', 'tau', 'beta', 'fat_R_B', 'v0', 'nx', 'R_ini', 'sigma_nominal', 'estado',
        'passo_falha', 'passos', 'sigma_med', 'umax_med', 'urms_med', 'umax_sigma', 'Ca_s', 'a4_med', 'a8_med',
        'K_med', 'neg_total_max', 'neg_cor_max', 'min_f', 'MLUPS', 'umax_rel_media', 'soma_Fx_med']

ORDEM_OPS = ['P1', 'P3 chi=-0.4', 'P3 chi=0.5', 'P3 chi=2', 'F1', 'F2', 'F3']


def num(x):
    try:
        return float(x)
    except ValueError:
        return float('nan')


def le_resumo(p):
    linhas = (p / 'resumo.dat').read_text().strip().splitlines()
    if len(linhas) < 2:
        return None
    v = linhas[1].split()
    d = dict(zip(COLS, v))
    for k in COLS:
        if k not in ('operador', 'estado'):
            d[k] = num(d.get(k, 'nan'))
    return d


def le_laplace(p):
    rows = []
    f = p / 'laplace.dat'
    if not f.exists():
        return rows
    for ln in f.read_text().splitlines():
        if ln.startswith('#') or not ln.strip():
            continue
        v = ln.split()
        rows.append([num(x) for x in v])
    return rows


def media(xs):
    xs = [x for x in xs if not math.isnan(x)]
    return sum(xs) / len(xs) if xs else float('nan')


def tendencia(rows, col=8):
    """col 8 ( zero-based ) = |u - v0|max no laplace.dat"""
    if len(rows) < 10:
        return float('nan')
    n = len(rows)
    fim = media([r[col] for r in rows[int(0.8 * n):]])
    meio = media([r[col] for r in rows[int(0.5 * n):int(0.7 * n)]])
    return fim / meio - 1.0 if meio > 0 else float('nan')


def coleta(raiz, lote):
    casos = []
    base = raiz / lote
    if not base.exists():
        return casos
    for r in sorted(base.rglob('resumo.dat')):
        p = r.parent
        d = le_resumo(p)
        if d is None:
            continue
        rows = le_laplace(p)
        d['pasta'] = str(p.relative_to(base))
        d['rotulo'] = d['pasta'].split('/')[0]
        d['tendencia'] = tendencia(rows)
        d['sigma_meiafinal'] = media([x[7] for x in rows[len(rows) // 2:]]) if rows else float('nan')
        d['R_meiafinal'] = media([x[1] for x in rows[len(rows) // 2:]]) if rows else float('nan')
        casos.append(d)
    return casos


def escreve_csv(casos, arq):
    if not casos:
        return
    chaves = ['pasta', 'rotulo'] + [k for k in casos[0] if k not in ('pasta', 'rotulo')]
    with open(arq, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=chaves, extrasaction='ignore')
        w.writeheader()
        for c in casos:
            w.writerow(c)


def e(x, d=3):
    return '—' if x is None or (isinstance(x, float) and math.isnan(x)) else f'{x:.{d}e}'


def ordena_ops(rotulos):
    return sorted(set(rotulos), key=lambda r: ORDEM_OPS.index(r) if r in ORDEM_OPS else 99)


#------------------------------------------------------------------------------------------------------------------#

def ajuste_quadratico(Rs, sig):
    """sigma = a + b / R + c / R^2  por minimos quadrados ( sem numpy ).  Devolve sigma_inf, delta, c2, R2."""
    X = [[1.0, 1.0 / R, 1.0 / R ** 2] for R in Rs]
    n = 3
    A = [[sum(X[k][i] * X[k][j] for k in range(len(X))) for j in range(n)] for i in range(n)]
    y = [sum(X[k][i] * sig[k] for k in range(len(X))) for i in range(n)]
    #  eliminacao de Gauss
    M = [A[i] + [y[i]] for i in range(n)]
    for i in range(n):
        piv = max(range(i, n), key=lambda r: abs(M[r][i]))
        M[i], M[piv] = M[piv], M[i]
        for r in range(i + 1, n):
            f = M[r][i] / M[i][i]
            for c in range(i, n + 1):
                M[r][c] -= f * M[i][c]
    b = [0.0] * n
    for i in reversed(range(n)):
        b[i] = (M[i][n] - sum(M[i][j] * b[j] for j in range(i + 1, n))) / M[i][i]
    a0, a1, a2 = b
    pred = [a0 + a1 / R + a2 / R ** 2 for R in Rs]
    m = media(sig)
    ss_res = sum((s - p) ** 2 for s, p in zip(sig, pred))
    ss_tot = sum((s - m) ** 2 for s in sig)
    return a0, a1 / a0, a2 / a0, 1 - ss_res / ss_tot if ss_tot > 0 else float('nan')


def tabela_E1(casos):
    out = ['## E1 — lei de Laplace ( τ = 1, σ nominal = 0.1 )', '',
           'σ( R ) = R Δp, média da metade final; ajuste σ∞ ( 1 + δ/R + c₂/R² ).', '',
           '| operador | R = 10 | 15 | 20 | 25 | 30 | 40 | 50 | σ∞ | δ | c₂ | R² |',
           '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
    for op in ordena_ops(c['rotulo'] for c in casos):
        cs = sorted([c for c in casos if c['rotulo'] == op and c['estado'] == 'estavel'], key=lambda c: c['R_ini'])
        por_R = {int(c['R_ini']): c for c in cs}
        cel = [f"{por_R[R]['sigma_meiafinal']:.6f}" if R in por_R else '—' for R in (10, 15, 20, 25, 30, 40, 50)]
        if len(cs) >= 4:
            s_inf, d, c2, r2 = ajuste_quadratico([c['R_meiafinal'] for c in cs], [c['sigma_meiafinal'] for c in cs])
            fit = [f'{s_inf:.6f}', f'{d:+.3f}', f'{c2:.3f}', f'{r2:.5f}']
        else:
            fit = ['—'] * 4
        out.append(f'| {op} | ' + ' | '.join(cel + fit) + ' |')
    return out


def tabela_E2(casos):
    out = ['## E2 — correntes espúrias ( R = 25 )', '',
           'Ca_s = ν |u|max / σ ( média da metade final ).  `*` = |u|max ainda variando mais de 5 %.  '
           '`X` = instável.', '']
    taus = sorted(set(c['tau'] for c in casos))
    sigs = sorted(set(c['sigma_nominal'] for c in casos))
    for op in ordena_ops(c['rotulo'] for c in casos):
        out += [f'### {op}', '', '| τ \\ σ | ' + ' | '.join(f'{s:g}' for s in sigs) + ' |',
                '|---|' + '---:|' * len(sigs)]
        for t in taus:
            cel = []
            for s in sigs:
                c = next((c for c in casos if c['rotulo'] == op and c['tau'] == t and c['sigma_nominal'] == s), None)
                if c is None:
                    cel.append('—')
                elif c['estado'] != 'estavel':
                    cel.append('X')
                else:
                    flag = '*' if abs(c['tendencia']) > 0.05 else ''
                    cel.append(e(c['Ca_s'], 2) + flag)
            out.append(f'| {t:g} | ' + ' | '.join(cel) + ' |')
        out.append('')
    return out


def tabela_E2chi(casos):
    out = ['## E2chi — varredura de χ ( P3, R = 25, σ = 0.1 )', '',
           '| τ | χ | \\|u\\|max/σ | a₄ | min f | tendência |', '|---:|---:|---:|---:|---:|---:|']
    for c in sorted(casos, key=lambda c: (c['tau'], c['chi'])):
        out.append(f"| {c['tau']:g} | {c['chi']:g} | {e(c['umax_sigma'])} | {e(c['a4_med'], 2)} | "
                   f"{c['min_f']:.4f} | {c['tendencia']:+.3f} |")
    return out


def tabela_E4(casos):
    out = ['## E4 — estabilidade da gota parada ( 30 000 passos )', '',
           '`ok` estável, `n` estável mas com populações totais negativas, `X` instável ( passo da falha ).', '']
    taus = sorted(set(c['tau'] for c in casos))
    sigs = sorted(set(c['sigma_nominal'] for c in casos))
    for op in ordena_ops(c['rotulo'] for c in casos):
        out += [f'### {op}', '', '| τ \\ σ | ' + ' | '.join(f'{s:g}' for s in sigs) + ' |',
                '|---|' + '---:|' * len(sigs)]
        for t in taus:
            cel = []
            for s in sigs:
                c = next((c for c in casos if c['rotulo'] == op and c['tau'] == t and c['sigma_nominal'] == s), None)
                if c is None:
                    cel.append('—')
                elif c['estado'] != 'estavel':
                    cel.append(f"X {int(c['passo_falha'])}")
                elif c['neg_total_max'] > 0:
                    cel.append('n')
                else:
                    cel.append('ok')
            out.append(f'| {t:g} | ' + ' | '.join(cel) + ' |')
        out.append('')
    return out


def tabela_E4v(casos):
    out = ['## E4v — gota em translação ( τ = 1, σ = 0.1 )', '',
           '\\|u − ⟨u⟩\\|max / σ ( relativa à velocidade média ), soma F_x por passo, populações de cor < 0.', '',
           '| operador | v0 | estado | \\|u−⟨u⟩\\|max/σ | Σ F_x | cor < 0 | total < 0 |', '|---|---:|---|---:|---:|---:|---:|']
    for op in ordena_ops(c['rotulo'] for c in casos):
        for c in sorted([c for c in casos if c['rotulo'] == op], key=lambda c: c['v0']):
            out.append(f"| {op} | {c['v0']:g} | {c['estado']} | {e(c['umax_rel_media'] / c['sigma_med'])} | "
                       f"{e(c['soma_Fx_med'], 2)} | {int(c['neg_cor_max'])} | {int(c['neg_total_max'])} |")
    return out


#------------------------------------------------------------------------------------------------------------------#
#   E3 -- interface inclinada ( resumo.dat de Imm_INCL, lido pelo cabecalho )
#------------------------------------------------------------------------------------------------------------------#

def coleta_incl(raiz, lote):
    casos = []
    base = raiz / lote
    if not base.exists():
        return casos
    for r in sorted(base.rglob('resumo.dat')):
        linhas = r.read_text().strip().splitlines()
        if len(linhas) < 2:
            continue
        cab = linhas[0].lstrip('#').split()
        d = {k: (v if k in ('operador', 'estado') else num(v)) for k, v in zip(cab, linhas[1].split())}
        d['pasta'] = str(r.parent.relative_to(base))
        d['rotulo'] = d['pasta'].split('/')[0]
        casos.append(d)
    return casos


def tabela_E3(casos):
    thetas = sorted(set(c['theta'] for c in casos))
    cab = '| operador | ' + ' | '.join(f'{t:.1f}°' for t in thetas) + ' |'
    sep = '|---|' + '---:|' * len(thetas)
    out = ['## E3 — interface plana inclinada ( N = 200, τ = 1, σ = 0.1 )', '']

    def bloco(titulo, f):
        linhas = [f'### {titulo}', '', cab, sep]
        for op in ordena_ops(c['rotulo'] for c in casos):
            cel = []
            for t in thetas:
                c = next((c for c in casos if c['rotulo'] == op and c['theta'] == t), None)
                cel.append('—' if c is None else ('X' if c['estado'] != 'estavel' else f(c)))
            linhas.append(f'| {op} | ' + ' | '.join(cel) + ' |')
        return linhas + ['']

    out += bloco('σ_pop / σ ( tensão nas populações; ~0 nas forças por construção )',
                 lambda c: f"{c['sigma_pop_med'] / c['sigma_nominal']:.5f}")
    out += bloco('espessura 10–90 % ( sítios )', lambda c: f"{2 * math.log(9) / c['K_med']:.4f}")
    out += bloco('\|u\|max / σ ( corrente espúria numa interface reta )', lambda c: e(c['umax/sigma_nominal'], 2))
    out += bloco('\|u_t\|max / \|u_n\|max', lambda c: f"{c['ut_max_med'] / c['un_max_med']:.2f}" if c['un_max_med'] > 0 else '—')
    out += bloco('deslocamento da interface d0 ( fim − início, sítios )', lambda c: f"{c['d0_fim'] - c['d0_ini']:+.1e}")
    return out


def tabela_E3chi(casos):
    out = ['## E3chi — interface inclinada θ = 26.6°, P3 contra χ', '',
           '| χ | σ_pop/σ | espessura | \|u\|max/σ | \|u_t\|/\|u_n\| |', '|---:|---:|---:|---:|---:|']
    for c in sorted(casos, key=lambda c: c['chi']):
        out.append(f"| {c['chi']:g} | {c['sigma_pop_med'] / c['sigma_nominal']:.5f} | "
                   f"{2 * math.log(9) / c['K_med']:.4f} | {e(c['umax/sigma_nominal'], 2)} | "
                   f"{c['ut_max_med'] / c['un_max_med'] if c['un_max_med'] > 0 else float('nan'):.2f} |")
    return out


TABELAS = {'E1': tabela_E1, 'E2': tabela_E2, 'E2chi': tabela_E2chi, 'E4': tabela_E4, 'E4v': tabela_E4v, 'E3': tabela_E3, 'E3chi': tabela_E3chi}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--saida', default='Lotes')
    a = ap.parse_args()
    raiz = (AQUI / a.saida) if not Path(a.saida).is_absolute() else Path(a.saida)

    md = ['# Resultados dos lotes', '']
    for lote, tab in TABELAS.items():
        casos = coleta_incl(raiz, lote) if lote.startswith('E3') else coleta(raiz, lote)
        if not casos:
            continue
        escreve_csv(casos, raiz / f'{lote}.csv')
        md += tab(casos) + ['']
        print(f'{lote}: {len(casos)} casos')

    (raiz / 'RESULTADOS.md').write_text('\n'.join(md) + '\n')
    print(f'escrito {raiz / "RESULTADOS.md"}')


if __name__ == '__main__':
    main()
