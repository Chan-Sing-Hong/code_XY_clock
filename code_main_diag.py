import argparse
from all_fn import *

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--model_name", type=str)
    parser.add_argument("--symm",  type=str)
    parser.add_argument("--T", type=float)
    parser.add_argument("-D", "--dcut",  type=int)
    parser.add_argument("-US", "--unit_size", type=int)
    parser.add_argument("--Lmax", type=int)
    parser.add_argument("--iL", type=int)
    parser.add_argument("--n", type=int)
    parser.add_argument("--num_collect", type=int)
    parser.add_argument("--hermitian", default=False, type=str)
    parser.add_argument("--verbose", default=False, type=str)

    args = parser.parse_args()

    model = Model(args.model_name, args.symm)
    T  = args.T
    
    HOTRG = Hotrg()
    HOTRG.dcut = args.dcut
    
    unit_size = args.unit_size
    Lmax = args.Lmax
    HOTRG.trgstep = int(np.floor(np.log2(Lmax/unit_size)))
    Ls = unit_size * 2**np.arange(0, HOTRG.trgstep+1)

    print("{}, symm={}, D={}".format(model.name, model.symm, HOTRG.dcut))
    print("T = ", T)
    print("L = ", Ls)
    
    L = Ls[args.iL]
    n = args.n
    num_collect = args.num_collect
    hermitian = args.hermitian

    print("D={}, n={}, num_collect={}".format(HOTRG.dcut,n,num_collect))

    # make folders
    prefix = f"{model.name}{'_symm' if model.symm else ''}/"
    Folder_Yi_raw = prefix + f"Yi_raw_data_wo_trans_num_collect{num_collect}/"
    Folder_Zi_raw = prefix + f"Zi_raw_data_wo_trans_num_collect{num_collect}/"
    makeFolder(Folder_Yi_raw)
    makeFolder(Folder_Zi_raw)

    T = args.T
    
    print('\nProcessing T = {}, L = {}'.format(T,L))
        
    Folder_T0 = prefix + f'T0_tensor/'
    T0 = cy.UniTensor.Load(Folder_T0 + "T0_T{}_L{}_D{}_US{}_{}.cytnx".format(T,L,HOTRG.dcut,unit_size,model.name))

    bi   = T0.bonds()[1]
    bi_n = bi.combineBonds([bi]*(n-1))
    d_x_n_sector_s = bi_n.getDegeneracies()

    if (n == 1) or (np.max(d_x_n_sector_s)-2 <= num_collect):
        print("Exact Diagonalization!")

        TM_n_traced = ExactCon.tTr_n([T0]*n)
        if hermitian:
            Y, Z = cy.linalg.Eigh(TM_n_traced)
        else:
            Y, Z = cy.linalg.Eig(TM_n_traced)
                
        Y.Save(Folder_Yi_raw + f"Y_T{T}_L{L}_D{HOTRG.dcut}_n{n}_US{unit_size}")
        Z.Save(Folder_Zi_raw + f"Z_T{T}_L{L}_D{HOTRG.dcut}_n{n}_US{unit_size}")

    else:
        for iq, qnum in enumerate(bi_n.qnums()):
            qnum = qnum[0]
            print("Partial Diagonalization!")
            d_x_n_sector = d_x_n_sector_s[iq]
            k  = np.min([num_collect, d_x_n_sector - 2])
            bo = cy.Bond(cy.BD_OUT, [[qnum]],[1], [cy.Symmetry.Zn(model.q)])

            phi_symm = cy.UniTensor([bi]*n + [bo], dtype=cy.Type.ComplexDouble)
            for isector, sector in enumerate(phi_symm.get_blocks()):
                phi_symm.put_block(cy.random.uniform(sector.shape(), low=-1., high=1., dtype=cy.Type.ComplexDouble), isector)
            
            LinOp = LinOp_n_PBC(d_x_n_sector, cy.Type.ComplexDouble, T0=T0, P0=None, qnum=qnum, n=n, get_momentum=False)
            
            try:            
                _ = cy.linalg.Lanczos(LinOp, k=k, Tin=phi_symm, which="LM", Maxiter=999999)
                Y = _[0]; Z_list = _[1:]

                Y.Save(Folder_Yi_raw + f"Y_qnum{qnum}_T{T}_L{L}_D{HOTRG.dcut}_n{n}_US{unit_size}")
                for iZ, Z in enumerate(Z_list):
                    Z.Save(Folder_Zi_raw + f"Z{iZ}_qnum{qnum}_T{T}_L{L}_D{HOTRG.dcut}_n{n}_US{unit_size}")
            except:
                print("partial diagonalization due to num_collect={} > d_x_n_sector={}".format(num_collect, d_x_n_sector))
                pass
    
    print("done")

if __name__ == "__main__":
    main()