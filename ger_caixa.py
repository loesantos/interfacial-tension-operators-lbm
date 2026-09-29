#!/usr/bin/env python3
"""Gera meio.vtk: caixa nx x ny x 1 toda fluida ( formato lido por read_geo ).

    python3 ger_caixa.py 100          # 100 x 100 x 1
    python3 ger_caixa.py 100 80       # 100 x 80 x 1
"""
import sys
nx = int(sys.argv[1]); ny = int(sys.argv[2]) if len(sys.argv) > 2 else nx
with open('meio.vtk', 'w') as f:
    f.write('# vtk DataFile Version 2.0\nGeometria\nASCII\nDATASET STRUCTURED_POINTS\n')
    f.write(f'DIMENSIONS {nx} {ny} 1\nASPECT_RATIO 1 1 1\nORIGIN 0 0 0\nPOINT_DATA {nx*ny}\n')
    f.write('SCALARS Geometria float\nLOOKUP_TABLE default\n')
    lin = ' '.join(['1'] * nx) + ' \n'
    for _ in range(ny): f.write(lin)
