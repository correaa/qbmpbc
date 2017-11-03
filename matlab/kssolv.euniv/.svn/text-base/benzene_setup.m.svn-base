%
% construct the Benzene molecule (ethane)
%

%
% 1. construct atoms
%
a1 = Atom('C');
a2 = Atom('H');
for j = 1:6
   atomlist(j) = a1;
end;
for j = 7:12
   atomlist(j) = a2;
end;
%atomlist = [a1 a1 a1 a1 a1 a1 ...
%            a2 a2 a2 a2 a2 a2];
%
% 2. set up the supercell
%
C = [20 0  0
     0 20  0
     0  0 10];
%
% 3. define the coordinates the atoms 
%
xyzlist = [
%
% coordinates of C atoms
%
  2.60782206   0.00000000   0.00000
  1.30391103   2.25844016   0.00000
 -1.30391103   2.25844016   0.00000
 -2.60782206   0.00000000   0.00000
 -1.30391103  -2.25844016   0.00000
  1.30391103  -2.25844016   0.00000
%
% coordinates of H atoms
%
  4.66762355   0.00000000   0.00000
  2.33381177   4.04228057   0.00000
 -2.33381177   4.04228057   0.00000
 -4.66762355   0.00000000   0.00000
 -2.33381177  -4.04228057   0.00000
  2.33381177  -4.04228057   0.00000

];
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',atomlist);
mol = set(mol,'xyzlist' ,xyzlist);
mol = set(mol,'ecut',25);  % kinetic energy cut off
mol = set(mol,'name','Benzene'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

