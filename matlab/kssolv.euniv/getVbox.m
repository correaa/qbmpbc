function [vBox] = getVbox(mol)
%
% Usage: [vBox] = getvbox(mol);
%
% Purpose:
%    Computes a box external potential.
%    V = 0 inside box
%      = 2000 meV outside
%    r0 is the center of the cell.
%
% Input:
%    mol       --- a Molecule object
%    L         --- size of box
% Output:
%    vBox --- potential (3D array)
%

global meanpotential

n1  = get(mol,'n1');
n2  = get(mol,'n2');
n3  = get(mol,'n3');

% Harmonic potential in real space.
vBox = 0/13.6056923*ones(n1,n2,n3);  %40eV
n1box = round(n1/2);
n2box = round(n2/2);
n3box = round(n3/2);
for i1 = 1:n1box
    for i2 = 1:n2box
        for i3 = 1:n3box
	    vBox(i1,i2,i3) = 0;
        end;
    end;
end;            

meanpotential=mean(mean(mean(vBox)));
disp(sprintf('mean potential = %e',meanpotential));

%vHarmonicf = fftn(vHarmonic);            

