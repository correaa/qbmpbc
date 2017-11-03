%
% construct the NiC molecule
%

%
% 1. construct atoms
%
a1 = Atom('Ni');
a2 = Atom('C');
atomlist(1) = a1;
atomlist(2) = a2;
%
% 2. set up the supercell
%
C = 5*eye(3);
%
% 3. define the coordinates the atoms in a.u.
%
xyzlist = [
0.0 0.0 0.0
1.8 0.0 0.0
]/0.529177;
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',atomlist);
mol = set(mol,'xyzlist' ,xyzlist);
mol = set(mol,'ecut',25);  % kinetic energy cut off
mol = set(mol,'name','NiC'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

