import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy as sp
import re
from pathlib import Path
import sys
# sys.path.insert(0, "/home/chansh/cytnx_install")
sys.path.insert(0, "/home/chansh/cytnx_install_12_Nov_2025")
import cytnx as cy

class Model():

    def __init__(self, name, symm):
        self.name = name
        self.symm = symm
        if name == "2D_Ising":
            self.q  = 2
            self.qnums = np.array([0,1])
            self.Tc = 2/np.log(1+np.sqrt(2))        
        elif re.match(r"2D_q=\d+_clock", name):
            q = int(re.search(r"q=(\d+)", name).group(1))         
            self.q  = q
            self.qnums = np.arange(0, q)
            self.Tc = {2: 2/np.log(np.sqrt(2)+1),
                       3: 3/(2*np.log(np.sqrt(3)+1)),
                       4: 1/np.log(np.sqrt(2)+1),
                       5: [0.90592, 0.95212],
                       6: [0.69014, 0.91275],
                       7: [0.53053, 0.90715],
                       8: [0.41723, 0.90605]}[q] # ref: arXiv:1912.11416
        elif re.match(r"2D_q=\d+_potts", name):
            q = int(re.search(r"q=(\d+)", name).group(1))
            self.q  = q
            self.qnums = np.arange(0, q)
            self.Tc = 1/np.log(1+np.sqrt(q))
        elif re.match(r"2D_M=\d+_XY", name):
            M = int(re.search(r"M=(\d+)", name).group(1))
            q = 2*M + 1
            qnums = np.arange(-M, M+1)
            qnums = np.where(qnums<0, qnums+q, qnums)
            self.M  = M            
            self.q  = q
            self.qnums = qnums
            self.Tc = 0.89301
        else:
            raise ValueError("Unknown model name")

    def λ_m(self, T, m):
        if self.name=="2D_Ising" or re.match(r"2D_q=\d+_clock", self.name):
            q = self.q        
            λ_m = 0
            for iq in range(q):
                θ_i = 2*np.pi / q * iq
                λ_m += np.cos(m * θ_i) * np.exp(1/T * np.cos(θ_i)) / np.sqrt(q)
            return λ_m
        elif re.match(r"(2D_q=\d+_potts)", self.name):
            q = self.q
            if m == 0:
                return (np.exp(1/T) + (q-1)) / np.sqrt(q)
            else:
                return (np.exp(1/T) - 1) / np.sqrt(q)
        elif re.match(r"2D_M=\d+_XY", self.name):
            M = self.M
            i = np.arange(-M, M+1)[m]
            # i = (m+1)//2 * (-1)**(m+1)
            return sp.special.iv(i, 1/T)

    def T0_1x1(self, T):
        q = self.q
        symm = self.symm
        T0_1x1 = cy.zeros((q,q,q,q))
        for i in range(q):
            for j in range(q):
                for k in range(q):
                    for l in range(q):
                        if (i+j-k-l) % q == 0:
                            T0_1x1[i,j,k,l] = np.sqrt(self.λ_m(T, i) * self.λ_m(T, j) * self.λ_m(T, k) * self.λ_m(T, l))
        
        if symm:
            bi = cy.Bond(cy.BD_IN,  [cy.Qs(qnum)>>1 for qnum in self.qnums], [cy.Symmetry.Zn(q)])
            bo = cy.Bond(cy.BD_OUT, [cy.Qs(qnum)>>1 for qnum in self.qnums], [cy.Symmetry.Zn(q)])
        else:
            bi  = cy.Bond(cy.BD_IN,  [cy.Qs(0)>>q], [cy.Symmetry.Zn(q)])
            bo  = cy.Bond(cy.BD_OUT, [cy.Qs(0)>>q], [cy.Symmetry.Zn(q)])
        T0_1x1_UT = cy.UniTensor([bi,bi,bo,bo])
        T0_1x1_UT.convert_from(cy.UniTensor(T0_1x1, rowrank=2))
        T0_1x1_UT.set_name("T0")
        T0_1x1_UT.set_labels(["u","l","r","d"])
        return T0_1x1_UT

    def P0_1x1(self):
            '''
            This function defines the OPE for 1-site translation operator P0
                                       0
                                  ┌────┘
                0                 │ ┏━━╳━━━┓ 
                │                 └─┨      ┃
            1───┐└────2    =    1───┨  P0  ┠──2
                │                   ┃      ┠─┐
                3                   ┗━━━━━━┛ │
                                    ┌─────┘
                                    3
            '''
            q = self.q
            symm = self.symm
            if symm == True:
                bi = cy.Bond(cy.BD_IN,  [cy.Qs(qnum)>>1 for qnum in self.qnums], [cy.Symmetry.Zn(q)])
                bo = cy.Bond(cy.BD_OUT, [cy.Qs(qnum)>>1 for qnum in self.qnums], [cy.Symmetry.Zn(q)])
            elif symm == False:
                bi = cy.Bond(cy.BD_IN,  [cy.Qs(0)>>q], [cy.Symmetry.Zn(q)])
                bo = cy.Bond(cy.BD_OUT, [cy.Qs(0)>>q], [cy.Symmetry.Zn(q)])
            else:
                raise TypeError("symm must be Boolean.")
            P0_UT = cy.UniTensor([bi,bi,bo,bo])

            P0 = cy.zeros((q,q,q,q))
            for i in range(q):
                for j in range(q):
                    for k in range(q):
                        for l in range(q):
                            if (i == k) & (j == l):
                                P0[i,j,k,l] = 1

            P0_UT.convert_from(cy.UniTensor(P0, rowrank=2))
            P0_UT.set_name("P0")
            P0_UT.set_labels(["u", "l", "r", "d"])

            return P0_UT

class ExactCon():

    def merge_two_x(Tl, Tr):
        Tl = Tl.clone().set_labels(["0","2","-1","4"])
        Tr = Tr.clone().set_labels(["1","-1","3","5"])
        Tl_Tr = cy.Contract(Tl, Tr)
        Tl_Tr.combineBonds(["0","1"], force=True)
        Tl_Tr.combineBonds(["4","5"], force=True)
        Tl_Tr.permute_(["0","2","3","4"])
        Tl_Tr.set_labels(["u", "l", "r", "d"])
        Tl_Tr.set_rowrank_(2)
        return Tl_Tr
    
    def merge_two_y(Tu, Td):
        Tu = Tu.clone().set_labels(["0","1","3","-1"])
        Td = Td.clone().set_labels(["-1","2","4","5"])
        Tu_Td = cy.Contract(Tu, Td)
        Tu_Td.combineBonds(["1","2"], force=True)
        Tu_Td.combineBonds(["3","4"], force=True)
        Tu_Td.permute_(["0","1","3","5"])
        Tu_Td.set_labels(["u", "l", "r", "d"])
        Tu_Td.set_rowrank_(2)
        return Tu_Td

    def contract_vertical(tensors):
        for i in range(len(tensors)):
            if i == 0:
                output = tensors[0]
            else:
                output = ExactCon.merge_two_y(output, tensors[i])
        return output
    
    def contract_horizontal(tensors):
        for i in range(len(tensors)):
            if i == 0:
                output = tensors[0]
            else:
                output = ExactCon.merge_two_x(output, tensors[i])
        return output

    def contract_LxL(tensors, unit_size):
        tensor_1xL_list = []
        for i in range(unit_size):
            tensors_vertical = [tensors[i+unit_size*j] for j in range(unit_size)]
            tensor_1xL = ExactCon.contract_vertical(tensors_vertical)
            tensor_1xL_list.append(tensor_1xL)

        tensor_LxL = ExactCon.contract_horizontal(tensor_1xL_list)
        return tensor_LxL
    
    def tTr_n(tensors):
        n = len(tensors)
        if n == 1:
            return tensors[0].Trace(0,3)
        else:
            for i,tensor in enumerate(tensors):
                if i == n-1:
                    tensor.set_labels([str(-(i+1)), str(2*i+1), str(2*i+101), str(-1)])
                else:
                    tensor.set_labels([str(-(i+1)), str(2*i+1), str(2*i+101), str(-(i+2))])

                # tensor.print_diagram()
                
                if i == 0:
                    output = tensor.clone()
                else:
                    output = cy.Contract(output, tensor)
                del tensor

            output.permute_([str(2*i+1) for i in range(n)] + [str(2*i+101) for i in range(n)])
            # output.combineBonds([str(2*i+1)   for i in range(n)], force=True)
            # output.combineBonds([str(2*i+101) for i in range(n)], force=True)
            output = output.set_rowrank(n)
        return output

    def matmul_n(Ml, Mr, n):
        Ml_Mr = cy.Contract(Ml.clone().set_labels([str(i)     for i in range(1,n+1)] + [str(100+i) for i in range(1,n+1)]),
                            Mr.clone().set_labels([str(100+i) for i in range(1,n+1)] + [str(200+i) for i in range(1,n+1)]))
        Ml_Mr.permute_([str(i) for i in range(1,n+1)] + [str(200+i) for i in range(1,n+1)])
        # Ml_Mr.combineBonds([str(i)     for i in range(1,n+1)], force=True)
        # Ml_Mr.combineBonds([str(200+i) for i in range(1,n+1)], force=True)
        Ml_Mr.set_rowrank_(n)
        return Ml_Mr

    def tTr_n_CBC(tensors):
        n = len(tensors)

        if n == 1:
            output = tensors[0].clone()
            output.set_labels(["2", "1", "101", "102"])
            output.permute_(["1", "2", "101", "102"])
            output.set_rowrank(2)
            return output
        else:
            for i in range(1,n+1):
                if i == 1:
                    labels = [str(n+1), str(i), str(100+i), str(-2*i)]
                elif i == n:
                    labels = [str(-2*(i-1)), str(i), str(100+i), str(101+n)]
                else:
                    labels = [str(-2*(i-1)), str(i), str(100+i), str(-2*i)]

                Ti = tensors[i-1].clone()
                Ti.set_labels(labels)
                if i == 1:
                    output = Ti
                else:
                    output = cy.Contract(output, Ti)
                    del Ti
            
            output.permute_([str(i) for i in range(1,n+2)] + [str(i+100) for i in range(1,n+2)])
            output.set_rowrank_(n+1)
            return output

class Hotrg():
    def __init__(self):
        self.name      = "HOTRG"
        self.trgstep   = None
        self.dcut      = None
        self.dcut_last = None

    def get_isometry(self, T1, T2, direction, dcut, min_blockdim=None, err=0):
        # get theta tensor
        if direction == "y":
            Tu = T1; Td = T2
            net = cy.Network("Networks/update_pure_y_symm_(cytnx_v1.0).net")
            net.PutUniTensor("Tu"  , Tu         , ["u", "l", "r", "d"])
            net.PutUniTensor("Tu.d", Tu.Dagger(), ["u", "l", "r", "d"])
            net.PutUniTensor("Td"  , Td         , ["u", "l", "r", "d"])
            net.PutUniTensor("Td.d", Td.Dagger(), ["u", "l", "r", "d"])
            net.PutUniTensors(["Tu", "Tu.d", "Td", "Td.d"],
                              [Tu, Tu.Dagger(), Td, Td.Dagger()])
        elif direction == "x":
            Tl = T1; Tr = T2
            net = cy.Network("Networks/update_pure_x_symm_(cytnx_v1.0).net")
            net.PutUniTensor("Tl"  , Tl         , ["u", "l", "r", "d"])
            net.PutUniTensor("Tl.d", Tl.Dagger(), ["u", "l", "r", "d"])
            net.PutUniTensor("Tr"  , Tr         , ["u", "l", "r", "d"])
            net.PutUniTensor("Tr.d", Tr.Dagger(), ["u", "l", "r", "d"])
            # net.PutUniTensors(["Tl", "Tl.d", "Tr", "Tr.d"],
            #                   [Tl, Tl.Dagger(), Tr, Tr.Dagger()])
        else:
            raise ValueError("direction must be x or y.")
        theta = net.Launch()
        
        # find isometry
        dc = np.min([dcut, np.prod(theta.shape()[:2]), np.prod(theta.shape()[2:])])
        if min_blockdim == None:
            S, U, Vd = cy.linalg.Svd_truncate(theta, keepdim=dc, err=err)
        else:
            S, U, Vd = cy.linalg.Svd_truncate(theta, keepdim=dc, min_blockdim=min_blockdim, err=err)
        del theta, Vd

        if direction == "y":
            U.set_labels(["x_u", "x_d", "dcut"])
        if direction == "x":
            U.set_labels(["y_l", "y_r", "dcut"])

        return S, U

    def update_pure(self, T1, T2, direction, dcut, Norm_divided, min_blockdim=None, err=0):
        # update tensors using isometry 
        if direction == "y":
            Tu = T1.clone(); Td = T2.clone()
            S, U = Hotrg.get_isometry(self, Tu, Td, direction, dcut, min_blockdim, err)

            net = cy.Network("Networks/update_two_y_symm_(cytnx_v1.0).net")
            net.PutUniTensor("Tu" , Tu        , ["u", "l", "r", "d"])
            net.PutUniTensor("Td" , Td        , ["u", "l", "r", "d"])
            net.PutUniTensor("U"  , U         , ["x_u", "x_d", "dcut"])
            net.PutUniTensor("U.d", U.Dagger(), ["x_u", "x_d", "dcut"])
            # net.PutUniTensors(['Tu', 'Td', 'U.d', 'U'], 
            #                   [Tu, Td, U.Dagger(), U])
            TT = net.Launch()
            del Tu, Td

        if direction == "x":
            Tl = T1.clone(); Tr = T2.clone()
            S, U = Hotrg.get_isometry(self, Tl, Tr, direction, dcut, min_blockdim, err)

            net = cy.Network("Networks/update_two_x_symm_(cytnx_v1.0).net")
            net.PutUniTensor("Tl" , Tl        , ["u", "l", "r", "d"])
            net.PutUniTensor("Tr" , Tr        , ["u", "l", "r", "d"])
            net.PutUniTensor("U"  , U         , ["y_l", "y_r", "dcut"])
            net.PutUniTensor("U.d", U.Dagger(), ["y_l", "y_r", "dcut"])
            # net.PutUniTensors(['Tl', 'Tr', 'U.d', 'U'],
            #                   [Tl, Tr, U.Dagger(), U])
            TT = net.Launch()
            del Tl, Tr

        TT.set_labels(["u", "l", "r", "d"])
        # normalization
        if Norm_divided == "Abs_Max":
            Norm = np.max([block.Abs().Max().item() for block in TT.get_blocks()])
        elif Norm_divided == "Norm":
            Norm = TT.Norm().item()
        elif Norm_divided == "Trace_T":
            Norm = TT.Trace("u", "d").Trace("l", "r").item()
        else:
            Norm = Norm_divided
        TT /= Norm

        return TT, U, Norm, S
    
    def put_ln_norm(self, Norms, N):
        ln_Ms = np.zeros_like(Norms)
        for i in range(N):
            ln_Ms[i] = np.log(Norms[i])*(4**(N-1-i))
        ln_norm_factor = np.sum(ln_Ms)
        return ln_norm_factor

    def update_impurity(self, T1, T2, direction, U, Norm_divided):
        # update tensors using isometry
        if direction == "y":
            Tu = T1.clone(); Td = T2.clone()
            net = cy.Network("Networks/update_two_y_symm_(cytnx_v1.0).net")
            net.PutUniTensor("Tu" , Tu        , ["u", "l", "r", "d"])
            net.PutUniTensor("Td" , Td        , ["u", "l", "r", "d"])
            net.PutUniTensor("U"  , U         , ["x_u", "x_d", "dcut"])
            net.PutUniTensor("U.d", U.Dagger(), ["x_u", "x_d", "dcut"])
            # net.PutUniTensors(['Tu', 'Td', 'U.d', 'U'], 
            #                   [Tu, Td, U.Dagger(), U])
            TT = net.Launch()
            del Tu, Td

        elif direction == "x":
            Tl = T1.clone(); Tr = T2.clone()
            net = cy.Network("Networks/update_two_x_symm_(cytnx_v1.0).net")
            net.PutUniTensor("Tl" , Tl        , ["u", "l", "r", "d"])
            net.PutUniTensor("Tr" , Tr        , ["u", "l", "r", "d"])
            net.PutUniTensor("U"  , U         , ["y_l", "y_r", "dcut"])
            net.PutUniTensor("U.d", U.Dagger(), ["y_l", "y_r", "dcut"])
            # net.PutUniTensors(['Tl', 'Tr', 'U.d', 'U'], 
            #                   [Tl, Tr, U.Dagger(), U])  
            TT = net.Launch()
            del Tl, Tr
        else:
            raise ValueError("direction must be x or y.")
    
        TT.set_labels(["u", "l", "r", "d"])
        # normalization
        if Norm_divided == "Abs_Max":
            Norm = np.max([block.Abs().Max().item() for block in TT.get_blocks()])
        elif Norm_divided == "Norm":
            Norm = TT.Norm().item()
        elif Norm_divided == "Trace_T":
            Norm = TT.Trace("u", "d").Trace("l", "r").item()
        else:
            Norm = Norm_divided
        TT /= Norm
    
        return TT, Norm

    def update_impurity_Tz(self, T1, T2, direction, U, Norm_divided):
        # update tensors using isometry
        if direction == "y":
            Tu = T1.clone(); Td = T2.clone()
            net = cy.Network("Networks/update_two_y_symm_Tz_(cytnx_v1.0).net")
            net.PutUniTensor("Tu" , Tu        , ["u", "l", "r", "d", "qnum"])
            net.PutUniTensor("Td" , Td        , ["u", "l", "r", "d", "qnum"])
            net.PutUniTensor("U"  , U         , ["x_u", "x_d", "dcut"])
            net.PutUniTensor("U.d", U.Dagger(), ["x_u", "x_d", "dcut"])
            # net.PutUniTensors(['Tu', 'Td', 'U.d', 'U'], 
            #                   [Tu, Td, U.Dagger(), U])
            TT = net.Launch()
            del Tu, Td

        elif direction == "x":
            Tl = T1.clone(); Tr = T2.clone()
            net = cy.Network("Networks/update_two_x_symm_Tz_(cytnx_v1.0).net")
            net.PutUniTensor("Tl" , Tl        , ["u", "l", "r", "d", "qnum"])
            net.PutUniTensor("Tr" , Tr        , ["u", "l", "r", "d", "qnum"])
            net.PutUniTensor("U"  , U         , ["y_l", "y_r", "dcut"])
            net.PutUniTensor("U.d", U.Dagger(), ["y_l", "y_r", "dcut"])
            # net.PutUniTensors(['Tl', 'Tr', 'U.d', 'U'], 
            #                   [Tl, Tr, U.Dagger(), U])  
            TT = net.Launch()
            del Tl, Tr
        else:
            raise ValueError("direction must be x or y.")

        TT.combineBonds(["qnum", "qnum_"], force=True)
        TT.set_labels(["u", "l", "r", "d", "qnum"])
        # normalization
        if Norm_divided == "Abs_Max":
            Norm = np.max([block.Abs().Max().item() for block in TT.get_blocks()])
        elif Norm_divided == "Norm":
            Norm = TT.Norm().item()
        elif Norm_divided == "Trace_T":
            Norm = TT.Trace("u", "d").Trace("l", "r").item()
        else:
            Norm = Norm_divided
        TT /= Norm
    
        return TT, Norm

    def run(self, T0, dcut, Norm_divided, get_momentum, P0):
        T0, Uy, Ny = Hotrg.update_pure(self, T0, T0, "y", dcut, Norm_divided)
        T0, Ux, Nx = Hotrg.update_pure(self, T0, T0, "x", dcut, Norm_divided)

        if get_momentum:
            P0, _ = Hotrg.update_impurity(self, P0, P0, "y", Uy, 1)

        return T0, Uy, Ux, Ny, Nx, P0

class Ptmrg(Hotrg):
    def __init__(self):
        super().__init__()
        self.name = "PTMRG"

    def run(self, T0_1x1, T0_1xL, T0_Lx1, T0_LxL, dcut, Norm_divided, get_momentum, P0_1x1, P0_1xL):

        T0_Lplus1xL, Ux, Norm_x1 = Ptmrg.update_pure    (self, T0_LxL, T0_1xL, "x", dcut, Norm_divided)
        T0_Lplus1x1,     Norm_x2 = Ptmrg.update_impurity(self, T0_Lx1, T0_1x1, "x", Ux  , Norm_divided)

        T0_Lplus1xLplus1, Uy, Norm_y1 = Ptmrg.update_pure    (self, T0_Lplus1xL, T0_Lplus1x1, "y", dcut, Norm_divided)
        T0_1xLplus1,          Norm_y2 = Ptmrg.update_impurity(self, T0_1xL     , T0_1x1     , "y", Uy  , Norm_divided)

        if get_momentum:
            P0_1xL, _ = Ptmrg.update_impurity(self, P0_1xL, P0_1x1, "y", Uy, 1)

        return T0_1xLplus1, T0_Lplus1x1, T0_Lplus1xLplus1, Uy, Ux, Norm_y1, Norm_y2, Norm_x1, Norm_x2, P0_1xL

class LinOp_n_PBC(cy.LinOp):
    def __init__(self, d_x_n_sector, dtype_cy, T0, P0, qnum, n, get_momentum):
        cy.LinOp.__init__(self, "mv", d_x_n_sector, dtype_cy)
        self.T0 = T0
        self.P0 = P0
        self.qnum = qnum
        self.n = n
        self.get_momentum = get_momentum

    def matvec(self, v0):
        output = v0
        output.set_labels([str(100+i) for i in range(1,self.n+1)] + [f"qnum{self.qnum}"])
        
        if self.get_momentum:  
            for i in range(1,self.n+1):
                Ti = self.P0.clone()
                if i != self.n:
                    Ti.set_labels([str(-2*i), str(i), str(i+100), str(-2*i-2)])
                else:
                    Ti.set_labels([str(-2*i), str(i), str(i+100), str(-2)])
                output = cy.Contract(output, Ti)
                del Ti
    
            output.permute_([str(i) for i in range(1,self.n+1)] + [f"qnum{self.qnum}"])
            output.set_labels([str(100+i) for i in range(1,self.n+1)] + [f"qnum{self.qnum}"])
            
        for i in range(1,self.n+1):
            Ti = self.T0.clone()
            if i != self.n:
                Ti.set_labels([str(-2*i), str(i), str(i+100), str(-2*i-2)])
            else:
                Ti.set_labels([str(-2*i), str(i), str(i+100), str(-2)])
            output = cy.Contract(output, Ti)
            del Ti

        output.permute_([str(i) for i in range(1,self.n+1)] + [f"qnum{self.qnum}"])
        return output

def makeFolder(Filename,sep='/'):
    '''file name may contains several folders, so need to check
    individually'''
    subFolderPath = Filename.split('/')
    complete = '.'
    for subFolderP in subFolderPath:
        complete = "/".join([complete,subFolderP])
        dir_path = Path(complete)
        if not dir_path.exists():
            dir_path.mkdir()

