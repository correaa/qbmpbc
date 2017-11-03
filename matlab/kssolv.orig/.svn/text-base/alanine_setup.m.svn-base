%
% construct the analine molecule
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
for j = 3:4
  atomlist(j) = a3;
end;
atomlist(5) = a4;
atomlist(6) = a3;
atomlist(7) = a1;
atomlist(8) = a4;
for j = 8:12
  atomlist(j) = a1;
end;

%atomlist = [a1; a2; a3; a3; a4; a3; a1; a4; a1; a1; a1; a1; a1];

%
% 2. set up the supercell
%
C = [20 0 0; 0 15 0; 0 0 20];
%
% 3. define the coordinates the atoms 
%
xyzlist = [
  -0.133895    1.636840    1.067415
  -0.114950    0.674034    1.390795
  -0.129352   -0.247512    0.272425
   1.116503   -0.074057   -0.581443
   1.816148    0.910827   -0.622102
  -1.393741   -0.013479   -0.566277
   0.745649    0.588173    1.923262
   1.378429   -1.162243   -1.308285
   2.122853   -1.195908   -1.879225
  -0.127349   -1.285733    0.685615
  -1.453821   -0.733080   -1.415204
  -1.409682    1.017453   -0.989385
  -2.313398   -0.144089    0.049763
]/0.529177;
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',atomlist);
mol = set(mol,'xyzlist' ,xyzlist);
mol = set(mol,'ecut',25);  % kinetic energy cut off
mol = set(mol,'name','alanine'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

