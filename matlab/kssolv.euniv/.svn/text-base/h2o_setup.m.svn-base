%
% construct the H2O molecule
%

%
% 1. construct atoms
%
a1 = Atom('O');
a2 = Atom('H');
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
% 3. define the coordinates the atoms in a.u.
%
coefs = [
0.000000000000    -0.123909374404     0.000000000000
1.429936611037     0.983265845431     0.000000000000
-1.429936611037     0.983265845431     0.000000000000
];
coefs = coefs/BoxSize; % Scale the atom coordinate relative to BoxSize.
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',atomlist);
mol = set(mol,'xyzlist' ,xyzlist);
mol = set(mol,'ecut',25);  % kinetic energy cut off
mol = set(mol,'name','H2O'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

