# acados benchmark
This repo contains my experiments with singularity free NMPC, as introduced in 
the publication: *NMPC for Racing Using a Singularity-Free Path-Parametric Model
 with Obstacle Avoidance - Daniel Kloeser, Tobias Schoels, Tommaso Sartor, Andrea Zanelli, Gianluca Frison, Moritz Diehl. Proceedings of the 21th IFAC World Congress, Berlin, Germany - July 2020*. 

This was used to benchmark acados on different machines using a relatively involved 
OCP formulation.

## workspace setup
Before running this make sure that you have installed cmake and python (e.g. using conda/mamba):
 
```bash
git clone https://github.com/acados/acados --filter=blob:none --recurse-submodules
cd acados
mkdir build && cd build
cmake -DACADOS_WITH_OPENMP=ON ..
make install
cd ..
export ACADOS_SOURCE_DIR=$(pwd)
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ACADOS_SOURCE_DIR/lib
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ACADOS_SOURCE_DIR/lib
pip3 install -e $ACADOS_SOURCE_DIR/interfaces/acados_template
pip3 install -U pip icecream
```
Now you can run the benchmark script using 

```bash
python3 main.py
```