%
% construct an electron quantum dot confined by 
% a harmonic external potential 
%
global ev meanpotential
%
% 2. set up the primitive cell 
%
C = 20*eye(3); % 40 Bohr = 2.1160 nm
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
vext = getV1overR(mol);
mol = set(mol,'vext', vext);
%mol = set(mol,'mpbc',1);
%
% set the number of electrons
%

mol = set(mol,'nel',0); 
mol = set(mol,'nspin',2);

X0 = genX0(mol,6);

options = setksopt;
options.X0 = X0;
% restart from a previous wave function
if (exist('X'))
    options.X0 = X;
end
options.maxcgiter = 100;
options.maxscfiter = 100;

%[Etotvec, X, vtot, rho] = scf(mol,options);
%plot_wf(X,1,mol);

% Representative results:
% L = 20; ecut = 20.5; mpbc = 0;
%   eigvals = [-0.4518   -0.1256   -0.1234   -0.1233   -0.1233   -0.0646];
% L = 20; ecut = 21; mpbc = 0;
%   eigvals = [-0.4540   -0.1258   -0.1234   -0.1233   -0.1233   -0.0647];
% L = 20; ecut = 21; mpbc = 1;
%   eigvals = [-0.4561   -0.1351   -0.1222   -0.1179   -0.1170   -0.0719];
% L = 15; ecut = 21; mpbc = 1;
%   eigvals = [-0.4544   -0.1507   -0.1141   -0.1103   -0.1091   -0.0569];
% L = 10; ecut = 21; mpbc = 1;
%   eigvals = [-0.4545   -0.1927   -0.0935   -0.0625   -0.0443    0.0093];
