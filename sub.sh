#!/bin/bash
#SBATCH --time=40:00:00   # walltime
#SBATCH --ntasks=1   # number of processor cores (i.e. tasks)
#SBATCH --nodes=1   # number of nodes
#SBATCH --cpus-per-task=2
 
#SBATCH --mem-per-cpu=3072M   # memory per CPU core
#SBATCH -J "dag"   # job name#
#./myscript.sh
  python mcts_cycle.py