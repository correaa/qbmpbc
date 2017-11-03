%
% construct the C12H26 molecule
%

%
% 1. construct atoms
%
a1 = Atom('C');
a2 = Atom('H');
atomlist(1) = a2;
atomlist(2) = a1;
for j = 3:5
   atomlist(j) = a2; 
end;  
for k = 0:10
   atomlist(5+3*k+1) = a1;
   atomlist(5+3*k+2) = a2;
   atomlist(5+3*k+3) = a2;
end;
%
%atomlist = [a2; a1; a2; a2; a2; a1; a2; a2; a1; a2; 
%            a2; a1; a2; a2; a1; a2; a2; a1; a2; a2; 
%            a1; a2; a2; a1; a2; a2; a1; a2; a2; a1;
%            a2; a2; a1; a2; a2; a1; a2; a2];
%
% 2. set up the supercell
%
C = [30 0 0; 0 15 0; 0 0 6];
%
% 3. define the coordinates the atoms 
%
xyzlist = [
   -5.4666    0.1256         0
   -5.6551   -0.9384         0
   -6.2646   -1.1498    0.8719
   -6.2646   -1.1498   -0.8719
   -4.4059   -2.4587    0.8599
   -4.3680   -1.7990         0
   -4.4059   -2.4587   -0.8599
   -2.4224   -1.4601    0.8648
   -2.9671   -1.1041         0
   -2.4224   -1.4601   -0.8648
   -3.5390    0.7544    0.8630
   -2.9671    0.4275         0
   -3.5390    0.7544   -0.8630
   -1.7159    1.9285    0.8560
   -1.6398    1.2669         0
   -1.7159    1.9285   -0.8560
    0.3226    1.1936    0.8630
   -0.1561    0.7490         0
    0.3226    1.1936   -0.8630
   -0.3226   -1.1936    0.8630
    0.1561   -0.7490         0
   -0.3226   -1.1936   -0.8630
    1.7159   -1.9285    0.8560
    1.6398   -1.2669         0
    1.7159   -1.9285   -0.8560
    3.5390   -0.7544    0.8630
    2.9671   -0.4275         0
    3.5390   -0.7544   -0.8630
    2.4224    1.4601    0.8648
    2.9671    1.1041         0
    2.4224    1.4601   -0.8648
    4.4059    2.4587    0.8599
    4.3680    1.7990         0
    4.4059    2.4587   -0.8599
    5.4666   -0.1256         0
    5.6551    0.9384         0
    6.2646    1.1498    0.8719
    6.2646    1.1498   -0.8719
]/0.529177;
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',atomlist);
mol = set(mol,'xyzlist' ,xyzlist);
mol = set(mol,'ecut',25);  % kinetic energy cut off
mol = set(mol,'name','c12h26'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

