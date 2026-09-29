#!/usr/bin/env python3
"""Roda os lotes de simulacao do artigo ( Interfacial Tension Paper ) em paralelo, na maquina local.

Cada caso e uma pasta com data_in.txt, operador.txt, velocity.txt e meio.vtk; o executavel e o
Imm_OPS desta pasta, rodado com UMA thread por caso e N casos ao mesmo tempo.  Um caso so conta como
feito quando tem resumo.dat: interrompido no meio, e refeito do zero na proxima chamada.  Assim o
lote pode ser parado ( Ctrl-C ) e retomado quando quiser.

Lotes
    E1      Lei de Laplace, R = 10 ... 50, caixa 4R, tau = 1, sigma = 0.1          ( calibracao )
    E2      correntes espurias, R = 25, tau x sigma                                ( E2 )
    E2chi   varredura fina de chi ( P3 ), R = 25, sigma = 0.1, tau = 0.6 e 1.0      ( E2 )
    E4      estabilidade da gota parada, tau x sigma alto, horizonte fixo           ( E4 )
    E4v     gota em translacao, v0 = 0.05 ... 0.25, tau = 1, sigma = 0.1            ( E4 )
    E3      interface plana inclinada ( Imm_INCL ), seis angulos, N = 200          ( E3 )
    E3chi   idem, varredura de chi ( P3 ) em theta = 26.6 graus ( p, q = 2, 1 )     ( E3 )
    todos   os sete, nessa ordem

Uso
    python3 rodar_lote.py E1 --procs 6
    python3 rodar_lote.py todos --procs 6 --lista        # so lista os casos e a estimativa de tempo
    python3 rodar_lote.py E2 --procs 6 --filtro F1        # so os casos cujo caminho contem 'F1'
    python3 rodar_lote.py E1 --procs 2 --passos 2000 --saida Lotes_teste      # teste rapido

Depois:  python3 analisar_lote.py  ( junta os resumos em tabelas )
"""

import argparse
import concurrent.futures as cf
import math
import os
import subprocess
import sys
import time
from pathlib import Path

AQUI = Path(__file__).resolve().parent
BIN = AQUI / 'Imm_OPS'
BIN_INCL = AQUI / 'Imm_INCL'

BETA = 0.8          # recoloracao de Latva-Kokko, largura 3.45 sitios
XI_F3 = 1.569988    # 2 / K em beta = 0.8
MLUPS_EST = 2.2     # uma thread, medido na maquina local ( teste/P3 e teste/F1v )

#  Operadores do nucleo da comparacao.  ( rotulo da pasta , linha do operador.txt , classe )

OPERADORES = [
    ('P1',            'P1',               'P'),
    ('P3 chi=-0.4',   'P3 -0.4',          'P'),     # Reis-Phillips
    ('P3 chi=0.5',    'P3 0.5',           'P'),     # Latva-Kokko & Rothman / Leclaire
    ('P3 chi=2',      'P3 2',             'P'),     # Liu-Valocchi-Kang
    ('F1',            'F1',               'F'),     # CSF
    ('F2',            'F2',               'F'),     # divergencia de tensao
    ('F3',            f'F3 {XI_F3}',      'F'),     # potencial quimico
]

CHI_FINO = [-2.0, -1.3, -1.0, -0.7, -0.55, -0.45, -0.4, -0.35, -0.3, -0.25, -0.15, 0.0, 0.5, 1.0, 2.0]


def nu(tau):
    return (tau - 0.5) / 3.0


def amplitude(classe, sigma, tau):
    """fat_R_B: A = 9 sigma / ( 4 tau ) nas perturbacoes, sigma nas forcas."""
    return 9.0 * sigma / (4.0 * tau) if classe == 'P' else sigma


def passos_convergencia(R, tau, minimo=30000, teto=100000):
    """Pelo menos duas vezes o tempo viscoso R^2 / nu, entre 30 000 e 100 000 passos."""
    return int(min(teto, max(minimo, math.ceil(2.0 * R * R / nu(tau) / 500.0) * 500)))


def caso(pasta, oper, classe, sigma, tau, nx, passos, v0=0.0, lote='', pq=None):
    return dict(pasta=pasta, oper=oper, classe=classe, sigma=sigma, tau=tau, nx=nx,
                passos=passos, v0=v0, lote=lote, A=amplitude(classe, sigma, tau), pq=pq,
                bin=BIN_INCL if pq else BIN)


def fmt(x):
    return f'{x:g}'


#------------------------------------------------------------------------------------------------------------------#
#   Definicao dos lotes
#------------------------------------------------------------------------------------------------------------------#

def lote_E1():
    casos = []
    for rot, oper, cl in OPERADORES:
        for R in (10, 15, 20, 25, 30, 40, 50):
            casos.append(caso(Path('E1') / rot / f'R={R}', oper, cl, 0.1, 1.0, 4 * R, 30000, lote='E1'))
    return casos


def lote_E2():
    casos = []
    for rot, oper, cl in OPERADORES:
        for tau in (0.6, 0.8, 1.0, 1.5):
            for sigma in (0.001, 0.003, 0.01, 0.03, 0.1, 0.3):
                casos.append(caso(Path('E2') / rot / f'tau={fmt(tau)}' / f'sigma={fmt(sigma)}',
                                  oper, cl, sigma, tau, 100, passos_convergencia(25, tau), lote='E2'))
    return casos


def lote_E2chi():
    casos = []
    for tau in (0.6, 1.0):
        for chi in CHI_FINO:
            casos.append(caso(Path('E2chi') / f'tau={fmt(tau)}' / f'chi={fmt(chi)}', f'P3 {chi}', 'P',
                              0.1, tau, 100, passos_convergencia(25, tau), lote='E2chi'))
    return casos


def lote_E4():
    #  Horizonte fixo de 30 000 passos: aqui a pergunta e "estoura ou nao", nao o valor convergido.
    casos = []
    for rot, oper, cl in OPERADORES:
        for tau in (0.51, 0.52, 0.55, 0.6, 0.8, 1.0, 1.5, 2.0):
            for sigma in (0.1, 0.3, 0.5, 1.0, 2.0):
                casos.append(caso(Path('E4') / rot / f'tau={fmt(tau)}' / f'sigma={fmt(sigma)}',
                                  oper, cl, sigma, tau, 100, 30000, lote='E4'))
    return casos


def lote_E4v():
    casos = []
    for rot, oper, cl in OPERADORES:
        for v0 in (0.05, 0.1, 0.15, 0.2, 0.25):
            casos.append(caso(Path('E4v') / rot / f'v0={fmt(v0)}', oper, cl, 0.1, 1.0, 100, 30000, v0=v0,
                              lote='E4v'))
    return casos


#  Interface plana inclinada: normal ( p, q ), caixa N x N.  Com N = 200 a menor distancia entre as duas
#  interfaces e 24 sitios ( p, q = 4, 1 ), sete espessuras.

INCLINACOES = [(1, 0), (4, 1), (3, 1), (2, 1), (3, 2), (1, 1)]     # 0, 14.0, 18.4, 26.6, 33.7, 45 graus


def lote_E3():
    casos = []
    for rot, oper, cl in OPERADORES:
        for p, q in INCLINACOES:
            casos.append(caso(Path('E3') / rot / f'pq={p},{q}', oper, cl, 0.1, 1.0, 200, 20000, lote='E3',
                              pq=(p, q)))
    return casos


def lote_E3chi():
    casos = []
    for chi in CHI_FINO:
        casos.append(caso(Path('E3chi') / f'chi={fmt(chi)}', f'P3 {chi}', 'P', 0.1, 1.0, 200, 20000,
                          lote='E3chi', pq=(2, 1)))
    return casos


LOTES = {'E1': lote_E1, 'E2': lote_E2, 'E2chi': lote_E2chi, 'E4': lote_E4, 'E4v': lote_E4v, 'E3': lote_E3, 'E3chi': lote_E3chi}


#------------------------------------------------------------------------------------------------------------------#
#   Montagem e execucao de um caso
#------------------------------------------------------------------------------------------------------------------#

DATA_IN = """meio.vtk

Dim_Pixel	1.0

steps:	    {passos}

files:	    0

n_threads:  1

tau:	    {tau}

tau_nh:     1.0

rho:	    1.0

tau_R:	    {tau}

tau_B:	    {tau}

tau_m:	    1.0

rho_R:	    1.0

rho_B:	    1.0

fat_R_B:	{A:.12g}

fat_R_R:	0.0

fat_B_B:    0.0

recoll:	    {beta}

wett_R:	    0.0
"""


def escreve_meio(pasta, nx):
    with open(pasta / 'meio.vtk', 'w') as f:
        f.write('# vtk DataFile Version 2.0\nGeometria\nASCII\nDATASET STRUCTURED_POINTS\n')
        f.write(f'DIMENSIONS {nx} {nx} 1\nASPECT_RATIO 1 1 1\nORIGIN 0 0 0\nPOINT_DATA {nx * nx}\n')
        f.write('SCALARS Geometria float\nLOOKUP_TABLE default\n')
        linha = ' '.join(['1'] * nx) + ' \n'
        f.write(linha * nx)


def monta(c, raiz):
    p = raiz / c['pasta']
    p.mkdir(parents=True, exist_ok=True)
    for velho in ('resumo.dat',):
        (p / velho).unlink(missing_ok=True)
    escreve_meio(p, c['nx'])
    (p / 'operador.txt').write_text(c['oper'] + '\n')
    (p / 'velocity.txt').write_text(fmt(c['v0']) + '\n')
    if c['pq']:
        (p / 'inclinacao.txt').write_text(f"{c['pq'][0]} {c['pq'][1]}\n")
    (p / 'data_in.txt').write_text(DATA_IN.format(passos=c['passos'], tau=c['tau'], A=c['A'], beta=BETA))
    return p


def feito(c, raiz):
    r = raiz / c['pasta'] / 'resumo.dat'
    return r.exists() and len(r.read_text().strip().splitlines()) >= 2


def roda(c, raiz):
    p = monta(c, raiz)
    t0 = time.time()
    with open(p / 'execucao.log', 'w') as log:
        ret = subprocess.run([str(c['bin'])], cwd=p, stdout=log, stderr=subprocess.STDOUT)
    estado = '?'
    if (p / 'resumo.dat').exists():
        linhas = (p / 'resumo.dat').read_text().strip().splitlines()
        if len(linhas) >= 2:
            cab = linhas[0].lstrip('#').split()
            v = linhas[1].split()
            estado = v[cab.index('estado')] if 'estado' in cab else "?"
    return c, ret.returncode, estado, time.time() - t0


def custo(c):
    return c['nx'] * c['nx'] * c['passos'] / (MLUPS_EST * 1e6)


def hms(s):
    s = int(s)
    return f'{s // 3600}h{(s % 3600) // 60:02d}m'


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('lote', choices=list(LOTES) + ['todos'])
    ap.add_argument('--procs', type=int, default=max(1, (os.cpu_count() or 2) // 2),
                    help='casos simultaneos ( padrao: metade dos nucleos logicos )')
    ap.add_argument('--saida', default='Lotes', help='pasta raiz dos casos ( padrao: Lotes )')
    ap.add_argument('--filtro', default='', help='so casos cujo caminho contem este texto')
    ap.add_argument('--passos', type=int, default=0, help='forca o numero de passos ( teste rapido )')
    ap.add_argument('--lista', action='store_true', help='so lista os casos e estima o tempo')
    ap.add_argument('--refaz', action='store_true', help='refaz tambem os casos ja concluidos')
    a = ap.parse_args()

    nomes = list(LOTES) if a.lote == 'todos' else [a.lote]
    casos = [c for n in nomes for c in LOTES[n]()]

    if not a.lista:
        for b in sorted(set(c['bin'] for c in casos)):
            if not b.exists():
                sys.exit(f'Nao achei {b}.  Compile antes:  make')
    if a.filtro:
        casos = [c for c in casos if a.filtro in str(c['pasta'])]
    if a.passos:
        for c in casos:
            c['passos'] = a.passos

    raiz = (AQUI / a.saida) if not Path(a.saida).is_absolute() else Path(a.saida)
    pend = [c for c in casos if a.refaz or not feito(c, raiz)]

    #  Os mais caros primeiro: equilibra melhor o fim do lote.
    pend.sort(key=custo, reverse=True)

    total = sum(custo(c) for c in pend)
    print(f'{len(casos)} casos no pedido, {len(pend)} a rodar, {a.procs} em paralelo')
    print(f'estimativa ( {MLUPS_EST} MLUPS por caso ): {hms(total)} de CPU, ~{hms(total / max(1, a.procs))} de relogio')
    print('( os instaveis param antes e custam menos )')

    if a.lista:
        for c in pend:
            print(f"   {str(c['pasta']):45s}  '{c['oper']}'  fat_R_B={c['A']:.6g}  passos={c['passos']}  ~{custo(c) / 60:.1f} min")
        return

    if not pend:
        print('nada a fazer.')
        return

    raiz.mkdir(parents=True, exist_ok=True)
    log = open(raiz / 'progresso.log', 'a')

    t0 = time.time()
    feitos, gasto = 0, 0.0
    ex = cf.ThreadPoolExecutor(max_workers=a.procs)
    try:
        futuros = [ex.submit(roda, c, raiz) for c in pend]
        for fu in cf.as_completed(futuros):
            c, rc, estado, dt = fu.result()
            feitos += 1
            gasto += custo(c)
            resta = (total - gasto) / max(1, a.procs)
            msg = (f'[{feitos}/{len(pend)}] {str(c["pasta"]):45s} {estado:9s} rc={rc}  {dt / 60:5.1f} min'
                   f'   decorrido {hms(time.time() - t0)}, resta ~{hms(resta)}')
            print(msg, flush=True)
            log.write(time.strftime('%Y-%m-%d %H:%M:%S ') + msg + '\n')
            log.flush()
        ex.shutdown(wait=True)
    except KeyboardInterrupt:
        #  O Ctrl-C do terminal chega tambem aos Imm_OPS em curso; os que estavam na fila sao cancelados.
        ex.shutdown(wait=False, cancel_futures=True)
        print('\nInterrompido.  Os casos concluidos ficam; rode de novo o mesmo comando para continuar.')
        sys.exit(1)


if __name__ == '__main__':
    main()
