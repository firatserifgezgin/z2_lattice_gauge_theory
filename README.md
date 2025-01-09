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

Once the representation of all operators are calculated, one can define a Hamiltonian and investigate the further properties of the system.

In the corresponding work, the Hamiltonian is,
$$H = -\cos\theta \sum_p U_p-\sin \theta \sum_i X_i.$$

Note that, $U_p=Z_i Z_j Z_k Z_l$ and the adiabatic parameter $\theta\in[0,\pi/2]$.

In the method part, we also provided dynamical correlation functions, entanglement entropy etc.
methods:

    % makespace_z2_lgt(n,m,mxSt): constructor, n x m plaq lattice, mxSt: maxString length 
    % o1.visStateM(idS{1},visualsOn,ifLabelPlaquettes) : visualize idS
    % [eigens,eigVecs] = o1.mbspectrum(theta)       : diag for all k-sec
    % [eigens,eigVecs] = o1.mbzerospectrum(theta)   : diag for k = 0
    % [eigens,eigVecs] = o1.mbspectrumKsector(theta,k) : diag for k = k0
    % [Xij]            = o1.S_XX_ij(obj,eigen,eigVec,i1,i2):S_XX(omega,i,j)
    % [Xkk_abs2,kk]    = o1.S_XX_kkp(eii,eiVec): S_XX(omega,k)
    % [Z4ij]           = o1.S_Z4Z4_ij(obj,eii,eiVec,i1,omega):S_UU(omega,i,j)
    % [Z4kk_abs,kk]    = o1.S_Z4Z4_kkp(eii,eiVec,i1,omega): S_UU(omega,k)
    % [S_A]            = o1.entEntropy(obj,theta,psi0) : calculate entanglement entropy for 

Clearly, we calculate the dynamical structure factors,
$$S_{XX}(\vec{k},\omega,\theta) = \langle X_{-\vec{k}}(\omega-H(\theta) + E_0 -i \eta)^{-1} X_{\vec{k}} \rangle_0$$
$$S_{UU}(\vec{k},\omega,\theta) = \langle U_{-\vec{k}}(\omega-H(\theta) + E_0 -i \eta)^{-1} U_{\vec{k}} \rangle_0$$
where $\eta$ is a small positive number which we use to smooth the spectra.
Here 
$$X_{\vec{k}}=\sum_j X_j e^{i{\vec k}\cdot \vec{r}_j},\quad
U_{\vec{k}}=\sum_p U_p e^{i{\vec k}\cdot \vec{R}_p}.$$
