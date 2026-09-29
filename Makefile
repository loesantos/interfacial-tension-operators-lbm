#=====================================================================================================================#
#   Makefile - interfacial-tension operators in colour-gradient LBM
#
#       make                builds Imm_OPS ( drop ), Imm_INCL ( inclined interface ) and verifica_operadores ( E0 )
#       make verifica       builds and runs the operator checks ( E0 ), writing E0_resultado.txt
#       make clean
#
#   SIMBOLTZ = folder with the SimBoltz functions ( default: ./SimBoltz_Functions, the subset included here )
#=====================================================================================================================#

SIMBOLTZ ?= SimBoltz_Functions
CXX      ?= g++
FLAGS    := -std=c++17 -O3 -march=native -fopenmp -I$(SIMBOLTZ)

DEPS     := Operadores_tensao.cpp Galilean_diagnostics.cpp

all: Imm_OPS Imm_INCL verifica_operadores

Imm_OPS: Imm_OPS.cpp $(DEPS)
	$(CXX) $(FLAGS) $< -o $@

Imm_INCL: Imm_INCL.cpp Operadores_tensao.cpp
	$(CXX) $(FLAGS) $< -o $@

verifica_operadores: verifica_operadores.cpp Operadores_tensao.cpp
	$(CXX) $(FLAGS) $< -o $@

verifica: verifica_operadores
	./verifica_operadores | tee E0_resultado.txt

clean:
	rm -f Imm_OPS Imm_INCL verifica_operadores

.PHONY: all verifica clean
