% 
%  construct bulk Silicon (crystal)
% 

% 
%  1. construct atoms
% 
a1 = Atom('Si');
atomlist = [a1; a1];
% 
%  2. set up the primitive cell 
% 
C = [
5.1306  5.1306   0.
5.1306   0.     5.1306
0.     5.1306   5.1306
] ; 
% 
%  3. define the coordinates the atoms
% 
coefs = [
0.0     0.0      0.0
0.25    0.25     0.25
] ; 
xyzlist = coefs*C';
% 
%  4. Configure the molecule (crystal)
% 
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',atomlist);
mol = set(mol,'xyzlist' ,xyzlist);
mol = set(mol,'ecut', 50);  % kinetic energy cut off
mol = set(mol,'name','Silicon');
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

