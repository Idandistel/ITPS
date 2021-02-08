import numpy as np
import scipy.sparse.linalg as sp
import itertools

def hyperelasticity(phase):

# ----------------------------------- GRID ------------------------------------

    ndim   = 3   # number of dimensions
    Nx      = 31  # number of voxels (assumed equal for all directions)
    Ny     =  1
    Nz=      1
    [Nx,Ny,Nz]=phase.shape
    shape  = [Nx,Ny,Nz]
    # ---------------------- PROJECTION, TENSORS, OPERATIONS ----------------------

    # tensor operations/products: np.einsum enables index notation, avoiding loops
    # e.g. ddot42 performs $C_ij = A_ijkl B_lk$ for the entire grid
    trans2 = lambda A2   : np.einsum('ijxyz          ->jixyz  ',A2   )
    ddot42 = lambda A4,B2: np.einsum('ijklxyz,lkxyz  ->ijxyz  ',A4,B2)
    ddot44 = lambda A4,B4: np.einsum('ijklxyz,lkmnxyz->ijmnxyz',A4,B4)
    dot11  = lambda A1,B1: np.einsum('ixyz   ,ixyz   ->xyz    ',A1,B1)
    dot22  = lambda A2,B2: np.einsum('ijxyz  ,jkxyz  ->ikxyz  ',A2,B2)
    dot24  = lambda A2,B4: np.einsum('ijxyz  ,jkmnxyz->ikmnxyz',A2,B4)
    dot42  = lambda A4,B2: np.einsum('ijklxyz,lmxyz  ->ijkmxyz',A4,B2)
    dyad22 = lambda A2,B2: np.einsum('ijxyz  ,klxyz  ->ijklxyz',A2,B2)
    print ("check")
    # identity tensor                                               [single tensor]
    i      = np.eye(ndim)
    # identity tensors                                            [grid of tensors]
    I      = np.einsum('ij,xyz'           ,                  i   ,np.ones([Nx,Ny,Nz]))
    I4     = np.einsum('ijkl,xyz->ijklxyz',np.einsum('il,jk',i,i),np.ones([Nx,Ny,Nz]))
    I4rt   = np.einsum('ijkl,xyz->ijklxyz',np.einsum('ik,jl',i,i),np.ones([Nx,Ny,Nz]))
    I4s    = (I4+I4rt)/2.
    II     = dyad22(I,I)

    # projection operator                                         [grid of tensors]
    # NB can be vectorized (faster, less readable), see: "elasto-plasticity.py"
    # - support function / look-up list / zero initialize


    Ghat4 = np.zeros([3, 3, 3, 3, Nx, Ny, Nz])  # projection operator
    x = np.zeros([3, Nx, Ny, Nz], dtype='int64')  # position vectors
    q = np.zeros([3, Nx, Ny, Nz], dtype='int64')  # frequency vectors
    delta = lambda i, j: np.float(i == j)  # Dirac delta function
    # - set "x" as position vector of all grid-points   [grid of vector-components]
    x[0], x[1], x[2] = np.mgrid[:Nx, :Ny, :Nz]
    # - convert positions "x" to frequencies "q"        [grid of vector-components]
    for i in range(3):
        freq = np.arange(-(shape[i] - 1) / 2, +(shape[i] + 1) / 2, dtype='int64')
        q[i] = freq[x[i]]
    # - compute "Q = ||q||", and "norm = 1/Q" being zero for the mean (Q==0)
    #   NB: avoid zero division
    q = q.astype(np.float)
    Q       = dot11(q,q)
    Z       = Q==0
    Q[Z]    = 1.
    norm    = 1./Q
    norm[Z] = 0.
    # - set projection operator                                   [grid of tensors]
    for i, j, l, m in itertools.product(range(3), repeat=4):
        Ghat4[i, j, l, m] = norm * delta(i, m) * q[j] * q[l]

        # (inverse) Fourier transform (for each tensor component in each direction)
        fft = lambda x: np.fft.fftshift(np.fft.fftn(np.fft.ifftshift(x), [Nx, Ny, Nz]))
        ifft = lambda x: np.fft.fftshift(np.fft.ifftn(np.fft.ifftshift(x), [Nx, Ny, Nz]))

        # functions for the projection 'G', and the product 'G : K^LT : (delta F)^T'
        G = lambda A2: np.real(ifft(ddot42(Ghat4, fft(A2)))).reshape(-1)
        K_dF = lambda dFm: trans2(ddot42(K4, trans2(dFm.reshape(3, 3, Nx, Ny, Nz))))
        G_K_dF = lambda dFm: G(K_dF(dFm))

        # ------------------- PROBLEM DEFINITION / CONSTITIVE MODEL -------------------

        # phase indicator: cubical inclusion of volume fraction (9**3)/(31**3)
        #phase  = np.zeros([N,N,N]); phase[-9:,:9,-9:] = 1.
        # material parameters + function to convert to grid of scalars
        param = lambda M0, M1: M0 * np.ones([Nx, Ny, Nz]) * (1. - phase) + \
                               M1 * np.ones([Nx, Ny, Nz]) * phase
        K      = param(0.833,8.33)  # bulk  modulus                   [grid of scalars]
        mu     = param(0.386,3.86)  # shear modulus                   [grid of scalars]

        # constitutive model: grid of "F" -> grid of "P", "K4"        [grid of tensors]
        def constitutive(F):
            C4 = K*II+2.*mu*(I4s-1./3.*II)
            S  = ddot42(C4,.5*(dot22(trans2(F),F)-I))            # Second piola kiroff stress
            P  = dot22(F,S)
            K4 = dot24(S,I4)+ddot44(ddot44(I4rt,dot42(dot24(F,C4),trans2(F))),I4rt)
            return P,K4

          # ----------------------------- NEWTON ITERATIONS -----------------------------

          # initialize deformation gradient, and stress/stiffness       [grid of tensors]
        F     = np.array(I,copy=True)
        P,K4  = constitutive(F)

        # set macroscopic loading
        DbarF = np.zeros([ndim,ndim,Nx,Ny,Nz]); DbarF[0,1] += 1.0

        # initial residual: distribute "barF" over grid using "K4"
        b     = -G_K_dF(DbarF)
        F    +=         DbarF
        Fn    = np.linalg.norm(F)
        iiter = 0

        # iterate as long as the iterative update does not vanish
        while True:
                dFm,_ = sp.cg(tol=1.e-8,
                  A = sp.LinearOperator(shape=(F.size,F.size),matvec=G_K_dF,dtype='float'),
                b = b,
                )                                        # solve linear system using CG
                F    += dFm.reshape(ndim,ndim,Nx,Ny,Nz)     # update DOFs (array -> tens.grid)
                P,K4  = constitutive(F)                  # new residual stress and tangent
                b     = -G(P)                            # convert res.stress to residual
                print('%10.2e'%(np.linalg.norm(dFm)/Fn)) # print residual to the screen
                if np.linalg.norm(dFm)/Fn<1.e-5 and iiter>0: break # check convergence
                iiter += 1

        def Stress_eq(F):
            C4 = K * II + 2. * mu * (I4s - 1. / 3. * II)
            S = ddot42(C4, .5 * (dot22(trans2(F), F) - I))
            s_mid=np.zeros([Nx,Ny,Nz])
            s_mid+=((np.square(S[0,0,:,:,:]-S[1,1,:,:,:])
            +np.square(S[1,1,:,:,:]-S[2,2,:,:,:])
            +np.square(S[2,2,:,:,:]-S[0,0,:,:,:])
            +6*(np.square(S[0,1,:,:,:])+np.square(S[1,2,:,:,:])+np.square(S[2,0,:,:,:])) / 2))
            print(s_mid.ndim)
            s_eq = np.sqrt(s_mid)                           # Von mises equivalent of S
            return s_eq

    s_eq = Stress_eq(F)
    print('finish calculation')
    return s_eq