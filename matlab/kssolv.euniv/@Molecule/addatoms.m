function mol = addatoms(mol, atomlist, xyzlist)
%
% usage: mol = add(mol, atomlist, xyzlist);
%
if (nargin == 3)
   %
   % check input
   %
   alist0       = mol.atomlist;
   mol.atomlist = [alist0; atomlist];
   xyzlist0     = mol.xyzlist;
   mol.xyzlist  = [xyzlist0; xyzlist];
   mol.natoms   = length(mol.atomlist);
else
   error('need exactly three arguements (mol, atomlist, position list)');
end;
