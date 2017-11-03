%
% construct the C2H6 molecule (ethane)
%

%
% 1. construct atoms
%
a1 = Atom('C');
a2 = Atom('H');
for j = 1:2
   atomlist(j) = a1;
end;
for j = 3:8
   atomlist(j) = a2;

end;
%atomlist = [a1; a1; a2; a2; a2; a2; a2; a2];
%
% 2. set up the supercell
%
BoxSize = 10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
   -0.000000000000    -0.000000000000    -1.417294491434
   -0.000000000000    -0.000000000000     1.417294491434
    1.603695532658    -0.000000000000    -1.983618366244
    0.801847766329     1.388841071218     1.983618366244
    0.801847766329    -1.388841071218     1.983618366244
   -1.603695532658     0.000000000000     1.983618366244
   -0.801847766329     1.388841071218    -1.983618366244
   -0.801847766329    -1.388841071218    -1.983618366244
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
mol = set(mol,'name','C2H6'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

