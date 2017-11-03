function [v1overR] = getV1overR(mol)
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

C = get(mol,'supercell');
r0 = 0.5*(C(:,1) + C(:,2) + C(:,3));
%r0 = r0 - [0.0 0.0 0.01]';

n1  = get(mol,'n1');
n2  = get(mol,'n2');
n3  = get(mol,'n3');
vol = get(mol,'vol');

% Harmonic potential in real space.
v1overR = [];
for i1 = 1:n1
    for i2 = 1:n2
        for i3 = 1:n3
            r = (i1-1)*C(:,1)/n1 + (i2-1)*C(:,2)/n2 + (i3-1)*C(:,3)/n3;            
            dr = r-r0;
            dr2 = norm(dr);
            if dr2~=0.0
                v1overR(i1,i2,i3) = -1/dr2;
            end
        end;
    end;
end;
v1overR(n1/2,n2/2,n3/2) = 2*v1overR(n1/2,n2/2,n3/2+1);
%v1overR = v1overR-max(max(max(v1overR)));
meanpotential=mean(mean(mean(v1overR)));
disp(sprintf('mean potential = %e',meanpotential));

%vHarmonicf = fftn(vHarmonic);            

