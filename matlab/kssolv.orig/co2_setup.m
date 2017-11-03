%
% construct the CO2 molecule
%

%
% 1. construct atoms
%
a1 = Atom('C');
a2 = Atom('O');
atomlist(1) = a1;
atomlist(2) = a2;
atomlist(3) = a2;
%atomlist = [a1; a2; a2];
%
% 2. set up the supercell
%
BoxSize = 10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
0.000000000000    -0.026381925226     0.000000000000
2.078381991574     0.009896367383     0.000000000000
-2.078381991574     0.009896367383     0.000000000000
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
mol = set(mol,'ecut', 25);  % kinetic energy cut off
mol = set(mol,'name','CO2'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

