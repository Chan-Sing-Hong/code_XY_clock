from all_fn import *

# Parameters setup
model_name = "#model_name#"
symm = #symm#
dcut = #dcut#
Lmax = #Lmax#
unit_size = #unit_size#
iteration = #iteration#
iT = #iT#

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
HOTRG.trgstep = int(np.floor(np.log2(Lmax/unit_size)))
Ls = unit_size * 2**np.arange(0, HOTRG.trgstep+1)

print("{}, symm={}, D={}".format(model.name, model.symm, HOTRG.dcut))
print("L = ", Ls)

# Folders creation
prefix = f"{model.name}{'_symm' if model.symm else ''}/"
Folder_T0    = prefix + f"T0_tensor/"
Folder_Norms = prefix + f"Norms/"
makeFolder(Folder_T0)
makeFolder(Folder_Norms)

for T in Ts:
    for L in Ls:
        print('\nProcessing T = {}, L = {}'.format(T,L))
        if L == unit_size:
            T0 = model.T0_1x1(T)
            T0 = ExactCon.contract_LxL([T0]*unit_size**2, unit_size)
            Norm = max([block.Abs().Max().item() for block in T0.get_blocks()])
            T0 /= Norm                
            Norms = [Norm]

        else:
            T0, Uy, Ny, Sy = HOTRG.update_pure(T0, T0, "y", HOTRG.dcut, "Abs_Max")
            del Uy, Sy
            T0, Ux, Nx, Sx = HOTRG.update_pure(T0, T0, "x", HOTRG.dcut, "Abs_Max")
            Norms = np.append(Norms, Ny**2 * Nx)
            del Ux, Sx, Ny, Nx

        T0.Save(Folder_T0 + "T0_T{}_L{}_D{}_US{}_{}".format(T,L,HOTRG.dcut,unit_size,model.name))
        np.save(Folder_Norms + "Norms_T{}_L{}_D{}_US{}_{}.npy".format(T,L,HOTRG.dcut,unit_size,model.name), Norms)

        print("done")