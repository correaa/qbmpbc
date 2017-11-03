%
% construct the SiH4 (Silane) molecule 
%

%
% 1. construct atoms
%
a1 = Atom('Si');
a2 = Atom('H');
alist(1) = a1;
for j = 2:5
  alist(j) = a2;
end; 
%alist = [a1; a2; a2; a2; a2];
%
% 2. set up supercell
%
C = 10*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [
 0.0     0.0      0.0
 0.161   0.161    0.161
-0.161  -0.161    0.161
 0.161  -0.161   -0.161
-0.161   0.161   -0.161
];
xyzmat = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',alist);
mol = set(mol,'xyzlist' ,xyzmat);
mol = set(mol,'ecut', 25);  % kinetic energy cut off
mol = set(mol,'name','SiH4'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;
