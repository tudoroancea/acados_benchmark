# acados benchmark
This repo contains my experiments with singularity free NMPC, as introduced in 
the publication: *NMPC for Racing Using a Singularity-Free Path-Parametric Model
 with Obstacle Avoidance - Daniel Kloeser, Tobias Schoels, Tommaso Sartor, Andrea Zanelli, Gianluca Frison, Moritz Diehl. Proceedings of the 21th IFAC World Congress, Berlin, Germany - July 2020*. 

This was used to benchmark acados on different machines using a relatively involved 
OCP formulation.

## workspace setup
Before running the benchmark make sure that you have installed 
- cmake and python, e.g. in a conda/mamba env using:
    ```bash 
    mamba create -n acados_benchmark python numpy matplotlib scipy casadi cmake compilers
    ```
- acados, e.g. using the following commands:
    ```bash
    git clone https://github.com/acados/acados --filter=blob:none --recurse-submodules
    cd acados
    mkdir build && cd build
    cmake -DACADOS_WITH_OPENMP=ON ..  # you can add other flags such as -DBLASFEO_TARGET=ARMV8A_ARM_CORTEX_A57
    make install
    cd ..
    export ACADOS_SOURCE_DIR=$(pwd)
    export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$ACADOS_SOURCE_DIR/lib
    export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:$ACADOS_SOURCE_DIR/lib
    pip3 install -e $ACADOS_SOURCE_DIR/interfaces/acados_template
    ```
- icecream (a python debugging logger), e.g. using 
  ```bash
  pip3 install icecream
  ```
  or 
  ```bash
  mamba install icecream
  ```

Now you can run the benchmark script using 

```bash
python3 main.py
```