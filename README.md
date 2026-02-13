# code_XY_clock
A research project implementing **HOTRG (Higher-Order Tensor Renormalization Group)** simulations using the [Cytnx](https://github.com/Cytnx-dev/Cytnx) library.

## 1. Environment Requirements
- **Operating System**: Linux or Windows Subsystem for Linux 2 (WSL2).

- **WSL2 Installation Guide**:\
[Official Microsoft WSL2 Instruction](https://learn.microsoft.com/en-us/windows/wsl/install/)

## 2. Anaconda Installation
You can download the installer from the [Anaconda Archive](https://repo.anaconda.com/archive/).

Or use the following command line:
```bash
curl -O https://repo.anaconda.com/archive/Anaconda3-2025.12-2-Linux-x86_64.sh
```

Install Anaconda:
```bash
bash ~/Anaconda3-2025.12-2-Linux-x86_64.sh
```

## 3. Cytnx Installation
Create a dedicated Conda environment build Cytnx from source.

### 3.1 Create Conda Environment
```bash
conda config --add channels conda-forge
conda create --name cytnx python=3.9 _openmp_mutex=*=*_llvm
conda activate cytnx
conda upgrade --all
```

### 3.2 Install Build Dependencies
```bash
conda install cmake make boost libboost git compilers numpy openblas pybind11 beartype boost-cpp arpack
```

### 3.3 Clone and Build
```bash
git clone https://github.com/cytnx-dev/cytnx.git
cd cytnx
mkdir build
cd build
```

Configure (Replace /home/chansh/cytnx_install with your desired path):
```bash
cmake .. -DCMAKE_INSTALL_PREFIX=/home/chansh/cytnx_install
```

Compile and install:
```bash
make -j$(nproc)
make install
```

## 4. Usage Configuration
**Path Setup**\
To use the compiled library in Python, add the installation path manually.

**Quick Start (In-script)**\
Add the following snippet at the beginning of your Python files (make sure this path is correctly modified in `all_fn.py`):
```python
import sys

# Update this path to where you installed Cytnx
sys.path.insert(0, "/home/chansh/cytnx_install")

import cytnx as cy
print(cy.__version__)
```

## 5. Code implement
### Parameters
`model_name`\
Defines the physical model
* 2D XY model, 
```python
model_name = "2D_M={}_XY"
```
Replace `{}` with the Fourier cutoff value $M$.

* 2D clock model:
```python
model_name = "2D_q={}_clock"
```
Replace `{}` with the q-value.

---
`symm`\
Boolean. Whether to preserve the symmetric block structure of tensors.

---

`dcut`\
The bond-dimension cutoff used in HOTRG.

---

`Lmax`\
Maximum system size 
$L$.

---

`unit_size`\
Controls the base system size.\
The initial $T$ tensor is exactly contracted to system size `unit_size * unit_size`.\
After N steps HOTRG, the system size becomes $L = $ unit_size $ \times 2^N $.

---

`iteration`\
The label of temperature list `Ts_tot`.

---

`iT`\
Index of temperature `T` within temperature list `Ts_tot`.

---

`iL`\
Index of size `L` within size list `Ls_tot`.

---

`n`\
Stacking number used to constructed transfer matrix TM $^{(L, nL)}$.\
The transfer matrix is built by tracing all vertical bonds using `n` copies of $T$ tensor.

Example: Transfer matrix (PBC) with `n=3`
````
          ┌───┐
        ┏━┷━┓ │
    1───┨ T ┠─│──101
        ┗━┯━┛ │                ┏━━╳━┓
        ┏━┷━┓ │              1─┨    ┠─101    
    2───┨ T ┠─│──102    =    2─┨ TM ┠─102
        ┗━┯━┛ │              3─┨    ┠─103
        ┏━┷━┓ │                ┗━━━━┛
    3───┨ T ┠─│──103
        ┗━┯━┛ │
          └───┘
````

---

`num_collect`\
Maximum number of eigenvalues to compute.

---

`hermitian`\
Boolean. Whether to perform Hermitian diagonalization.

---