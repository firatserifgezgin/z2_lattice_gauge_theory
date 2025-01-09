# z2_lattice_gauge_theory
We keep our numerical code repository pertaining to your manuscript arxiv: 2408.14295

It is a MATLAB implementation of a Z2 lattice gauge theory. 

The main file is the class desccription of a 'Z2 object'. The properties of this object could be accesses as:

properties:
    % o1.idS        : generates gauge-invariant basis (all loops of size X)
    % o1.xmat,zmat  : representation X and Z in idS
    % o1.X_k        : representation of X_k for each k-sector
    % o1.Z4         : representation of ZZZZ for each plaquette
    % o1.ksec       : clustering of idS into k-sectors and projectors
    % methods:
    % makespace_z2_lgt(n,m,mxSt): constructor, n x m plaq lattice,
    %                                       mxSt: maxString length 
    % o1.visStateM(idS{1},visualsOn,ifLabelPlaquettes) : visualize idS
    % [eigens,eigVecs] = o1.mbspectrum(theta)       : diag for all k-sec
    % [eigens,eigVecs] = o1.mbzerospectrum(theta)   : diag for k = 0
    % [eigens,eigVecs] = o1.mbspectrumKsector(theta,k) : diag for k = k0
    % [Xij]            = o1.S_XX_ij(obj,eigen,eigVec,i1,i2):S_XX(omega,i,j)
    % [Xkk_abs2,kk]    = o1.S_XX_kkp(eii,eiVec): S_XX(omega,k)
    % [Z4ij]           = o1.S_Z4Z4_ij(obj,eii,eiVec,i1,omega):S_UU(omega,i,j)
    % [Z4kk_abs,kk]    = o1.S_Z4Z4_kkp(eii,eiVec,i1,omega): S_UU(omega,k)


