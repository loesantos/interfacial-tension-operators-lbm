#!/bin/bash
# roda_caso.sh <pasta> "<operador>" <A ou sigma> <tau> <beta> <passos> <nx> [v0]
# ex.: ./roda_caso.sh teste/P3 "P3 0.5" 0.225 1.0 0.8 30000 100 0
# Na classe P, sigma = (4/9) A tau; na classe F, o terceiro argumento ja e sigma.
set -e
dir="$1"; oper="$2"; amp="$3"; tau="$4"; beta="$5"; steps="$6"; nx="$7"; v0="${8:-0}"
BIN="$(cd "$(dirname "$0")" && pwd)/Imm_OPS"
GER="$(cd "$(dirname "$0")" && pwd)/ger_caixa.py"
mkdir -p "$dir"; cd "$dir"
python3 "$GER" "$nx"
echo "$oper" > operador.txt
echo "$v0" > velocity.txt
cat > data_in.txt << EOD
meio.vtk

Dim_Pixel	1.0

steps:	    $steps

files:	    0

n_threads:  1

tau:	    $tau

tau_nh:     1.0

rho:	    1.0

tau_R:	    $tau

tau_B:	    $tau

tau_m:	    1.0

rho_R:	    1.0

rho_B:	    1.0

fat_R_B:	$amp

fat_R_R:	0.0

fat_B_B:    0.0

recoll:	    $beta

wett_R:	    0.0
EOD
"$BIN" > execucao.log 2>&1 || true
tail -1 resumo.dat
