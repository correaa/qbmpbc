%
% construct an electron quantum dot confined by 
% a harmonic external potential 
%

%
% 2. set up the primitive cell 
%
 C = [10  0  0;
       0 10  0;
       0  0  0.5]; % 40 Bohr = 2.1160 nm
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'ecut', 21);  % kinetic energy cut off
mol = set(mol,'name','Quantum dot'); 
n1 = get(mol,'n1');
n2 = get(mol,'n2');
n3 = get(mol,'n3');
% 
% construct external potential
%
vext = getVbox(mol); % zero potential
mol = set(mol,'vext', vext);
mol = set(mol,'mpbc',1);
%
% set the number of electrons
%

mol = set(mol,'nel',0); 
mol = set(mol,'nspin',2);

X0 = genX0(mol,6);

options = setksopt;
options.X0 = X0;
options.maxcgiter  = 100;
options.maxscfiter = 100;

global ev

%[Etotvec, X, vtot, rho] = scf(mol,options);
%plot_wf(X,1,mol);
