function totnel = getnel(mol)
%
% count the total number of valence electrons in the molecule (mol)
%
% usage: totnel = getnel(mol);
%

atomlist = get(mol,'atomlist');
if (~isempty(atomlist))
   natoms = get(mol,'natoms');
   totnel = 0;
   for j = 1:natoms
      atom = atomlist(j);
      totnel = totnel + get(atom,'venum');
   end;
else
   totnel = get(mol,'nel');
end;
