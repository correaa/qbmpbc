%
% construct an Aluminum cluster
%

%
% 1. construct atoms
%
a = Atom('Al');
for j =1:8
   alist(j) = a;
end;
%
% 2. set up supercell
%
C = 20*eye(3); 
%
% 3. define the coordinates the atoms 
%

xyzmat = [
0 0 0
10 0 0
0 10 0
0 0 10
10 10 10
10 10 0
0 10 10
10 0 10
];
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',alist);
mol = set(mol,'xyzlist' ,xyzmat);
mol = set(mol,'ecut', 25);  % kinetic energy cut off
mol = set(mol,'name','Aluminum cluster'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;
