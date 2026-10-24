#!/bin/bash
set -euo pipefail
project_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)
source "$project_dir/scripts/activate.sh"
mkdir -p "$project_dir/runs"
run=$(mktemp -d "$project_dir/runs/water-smoke.XXXXXX")
cd "$run"
cat > topol.top <<'TOP'
#include "amber99sb-ildn.ff/forcefield.itp"
#include "amber99sb-ildn.ff/tip3p.itp"
[ system ]
Water environment smoke test
[ molecules ]
TOP
gmx solvate -cs spc216.gro -box 2 2 2 -o water.gro -p topol.top > solvate.log 2>&1
cat > em.mdp <<'MDP'
integrator = steep
nsteps = 100
emtol = 1000
cutoff-scheme = Verlet
rlist = 0.8
coulombtype = PME
rcoulomb = 0.8
rvdw = 0.8
pbc = xyz
MDP
gmx grompp -f em.mdp -c water.gro -p topol.top -o em.tpr > em-prep.log 2>&1
gmx mdrun -deffnm em -ntmpi 1 -ntomp 2 > em-run.log 2>&1
cat > md.mdp <<'MDP'
integrator = md
nsteps = 100
dt = 0.002
cutoff-scheme = Verlet
rlist = 0.8
coulombtype = PME
rcoulomb = 0.8
rvdw = 0.8
constraints = h-bonds
tcoupl = v-rescale
tc-grps = System
tau-t = 0.1
ref-t = 300
pcoupl = no
gen-vel = yes
gen-temp = 300
gen-seed = 42
nstenergy = 10
pbc = xyz
MDP
gmx grompp -f md.mdp -c em.gro -p topol.top -o md.tpr > md-prep.log 2>&1
gmx mdrun -deffnm md -ntmpi 1 -ntomp 2 > md-run.log 2>&1
gmx check -e md.edr > check.log 2>&1
test -s md.gro
printf 'Water smoke test passed: %s\n' "$run"
