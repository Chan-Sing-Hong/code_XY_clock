from all_fn import *

model_name = "#model_name#"
symm = #symm#
dcut = #dcut#
Lmax = #Lmax#
unit_size = #unit_size#
iteration = #iteration#
iT = #iT#
iL = #iL#
n = #n#
num_collect = #num_collect#
hermitian = #hermitian#

model = Model(model_name, symm)

if iteration == 300:
    Ts_tot = model.Tc
    if type(Ts_tot) != list:
        Ts_tot = [Ts_tot]
elif iteration == 0:
    Ts_tot = np.arange(0.8, 0.96+0.02, 0.02)

Ts = [Ts_tot[iT]]
# Ts = Ts_tot

HOTRG = Hotrg()
HOTRG.dcut = dcut
trgstep = int(np.floor(np.log2(Lmax/unit_size)))
Ls_tot = unit_size * 2**np.arange(0, trgstep+1)
Ls = [Ls_tot[iL]]
# Ls = Ls_tot

print("{}, symm={}, D={}".format(model.name, model.symm, HOTRG.dcut))
print("L = ", Ls)
print("D={}, n={}, num_collect={}".format(HOTRG.dcut,n,num_collect))

# make folders
prefix = f"{model.name}{'_symm' if model.symm else ''}/"
Folder_Yi_raw = prefix + f"Yi_raw_data_wo_trans_num_collect{num_collect}/"
Folder_Zi_raw = prefix + f"Zi_raw_data_wo_trans_num_collect{num_collect}/"
makeFolder(Folder_Yi_raw)
makeFolder(Folder_Zi_raw)

for T in Ts:
    for L in Ls:
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

            print("done")

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

                    print("done")
                except:
                    print("partial diagonalization due to num_collect={} > d_x_n_sector={}".format(num_collect, d_x_n_sector))
                    pass