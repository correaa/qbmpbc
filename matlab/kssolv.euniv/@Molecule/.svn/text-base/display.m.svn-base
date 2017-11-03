function display(mol)
%
% display method for atom
%
fprintf('Molecule: %s\n', mol.name);
fprintf('   supercell:\n');
if (~isempty(mol.supercell))
   fprintf('    %11.3e   %11.3e   %11.3e\n', ...
           mol.supercell(1,1), mol.supercell(1,2), mol.supercell(1,3));
   fprintf('    %11.3e   %11.3e   %11.3e\n', ...
           mol.supercell(2,1), mol.supercell(2,2), mol.supercell(2,3));
   fprintf('    %11.3e   %11.3e   %11.3e\n', ...
           mol.supercell(3,1), mol.supercell(3,2), mol.supercell(3,3));
end;
fprintf('   sampling size: n1 = %d, n2 = %d, n3 = %d\n', mol.n1, mol.n2, mol.n3);
fprintf('   atoms and coordinates: ');
alist   = mol.atomlist;
if (~isempty(alist))
    xyzlist = mol.xyzlist;
    fprintf('\n');
    for i = 1:mol.natoms
       a = alist(i);
       xyz = xyzlist(i,:);
       symbol = get(a,'symbol');
       fprintf('   %5d   %2s    %11.3e    %11.3e    %11.3e\n', ...
       i, symbol, xyz(1), xyz(2), xyz(3));
    end; 
else
    fprintf('none\n'); 
end;
fprintf('   number of electrons  : %d\n', mol.nel);
fprintf('   spin type            : %d\n', mol.nspin);
fprintf('   kinetic energy cutoff: %9.3e\n', mol.ecut);
if (~isempty(alist))
   fprintf('   pseudo pot cutoff    : %9.3e (A.U.)\n', mol.rcut);
else
   fprintf('   pseudo pot cutoff    : none \n');
end;
fprintf('   temperature: %9.3e\n', mol.temperature);
if (mol.mpbc ~= 0)
    fprintf('   mpbc: %d, mpbc_denom: %d\n',mol.mpbc,mol.mpbc_denom);
end
