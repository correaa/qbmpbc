%
% construct the planar singlet silysilylene Si2H4 molecule (crystal)
%

%
% 1. construct atoms
%  J. Am. Chem. Soc.  p5222 (1979)

a1 = Atom('Si');
a2 = Atom('H');
%atomlist = [a1; a1; a2; a2; a2; a2];
for j = 1:2
  atomlist(j) = a1;
end;
for j = 3:6
  atomlist(j) = a2;
end;

%atomlist = [a1; a1; a2; a2; a2; a2];
len_SiSi = 2.083;      % Si-Si bond length
len_SiH = 1.482;       % Si-H bond length
angle = 114.4/180*pi;  % H-Si-H angle

%
% 2. set up the supercell
%
BoxSize=10;
C = BoxSize*eye(3);
%
% 3. define the coordinates the atoms 
%
coefs = [      
-len_SiSi*0.5  0.0   0.0
 len_SiSi*0.5  0.0   0.0
 len_SiSi*0.5+len_SiH*cos(angle/2)  len_SiH*sin(angle/2)  0.0  
 len_SiSi*0.5+len_SiH*cos(angle/2) -len_SiH*sin(angle/2)  0.0
-len_SiSi*0.5-len_SiH*cos(angle/2)  len_SiH*sin(angle/2)  0.0  
-len_SiSi*0.5-len_SiH*cos(angle/2) -len_SiH*sin(angle/2)  0.0  
];
coefs = coefs/BoxSize;
xyzlist = coefs*C';
%
% 4. Configure the molecule (crystal)
%
mol = Molecule();
mol = set(mol,'supercell',C);
mol = set(mol,'atomlist',atomlist);
mol = set(mol,'xyzlist' ,xyzlist);
mol = set(mol,'ecut', 25);  % kinetic energy cut off
mol = set(mol,'name','Si2H4'); 
ierr = finalize(mol);
if (ierr ~=0)
   fprintf('Error in the molecule setup\n');
end;

