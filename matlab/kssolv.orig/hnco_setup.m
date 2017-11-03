%
% construct the HNCO (isocyanic acid) molecule
%

%
% 1. construct atoms
%
a1 = Atom('H');
a2 = Atom('N');
a3 = Atom('C');
a4 = Atom('O');
atomlist(1) = a1;
atomlist(2) = a2;
atomlist(3) = a3;
atomlist(4) = a4;

%atomlist = [a1; a2; a3; a4];
%
% 2. set up the supercell
%
BoxSize = 10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
 3.284938618980    -1.461359490637     0.001371230530
 2.258342591219     0.125197975455     0.001932017162
-0.006404171175     0.010271400781    -0.004893948816
-2.179288390874    -0.025234185924     0.001893804798
];
coefs = coefs/BoxSize;
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',atomlist);
mol = set(mol,'xyzlist' ,xyzlist);
mol = set(mol,'ecut',25);  % kinetic energy cut off
mol = set(mol,'name','HNCO'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

