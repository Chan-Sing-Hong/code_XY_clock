# code_XY_clock
A research project utilizing the Cytnx library for HOTRG simulations.

## 1. Environment Requirements
- **Operating System**: Linux or Windows Subsystem for Linux 2 (WSL2).

- **WSL2 Installation Guide**: [Official Microsoft WSL2 Instruction](https://learn.microsoft.com/en-us/windows/wsl/install/)

## 2. Anaconda Installation
You can download the installer from the [Anaconda Archive](https://repo.anaconda.com/archive/) or use the following command-line steps:
````
curl -O https://repo.anaconda.com/archive/Anaconda3-2025.12-2-Linux-x86_64.sh
````

Install it with
````
bash ~/Anaconda3-2025.12-2-Linux-x86_64.sh
````

## 3. Cytnx Installation
Follow these steps to set up a dedicated Conda environment and compile Cytnx.

### 3.1 Create Conda Environment
````
conda config --add channels conda-forge
conda create --name cytnx python=3.9 _openmp_mutex=*=*_llvm
conda activate cytnx
conda upgrade --all
````

### 3.2 Install Build Dependencies
````
conda install cmake make boost libboost git compilers numpy openblas pybind11 beartype boost-cpp arpack
````

### 3.3 Clone and Build
````
git clone https://github.com/cytnx-dev/cytnx.git
cd cytnx
mkdir build
cd build

# Configure (Replace /home/chansh/cytnx_install with your desired path)
cmake .. -DCMAKE_INSTALL_PREFIX=/home/chansh/cytnx_install

# Compile and Install
make -j$(nproc)
make install
````

## 4. Usage Configuration
**Path Setup**\
To use the compiled library in your Python scripts, you must point to the installation directory.

**In-script (Quick Start)**\
Add the following snippet at the beginning of your Python files (ensure this is modified in all_fn.py):
````
import sys

# Update this path to where you installed Cytnx
sys.path.insert(0, "/home/chansh/cytnx_install")

import cytnx as cy
print(cy.__version__)
````
Please modify this path in `all_fn.py`.

