%
% construct an electron quantum dot confined by 
% a harmonic external potential 
%

%
% 2. set up the primitive cell 
%
C = 10*eye(3);
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'ecut', 25);  % kinetic energy cut off
mol = set(mol,'name','Quantum dot'); 
% 
% construct external potential
%
vext = getVharmonic(mol,1);
mol  = set(mol,'vext', vext);
%
% set the number of electrons
%
mol = set(mol,'nel',8); 
mol = set(mol,'nspin',2);

