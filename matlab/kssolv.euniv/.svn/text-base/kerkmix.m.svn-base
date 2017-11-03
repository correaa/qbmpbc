function vnew = kerkmix(mol, vin, vout)
%
% Usage: vnew = kerkmix(mol, vin, vout);
%
% Purpose: 
%     Construct a new 3D potential by mixing the potentials vin and vout
%     using Kerker's recipe;
%
% Input:
%     mol  - a Molecule object
%     vin  - a 3D potential
%     vout - a 3D potential
%
% Output:
%     vnew - a 3D output potential  
%
[n1,n2,n3]=size(vin);
n123 = n1*n2*n3;
vres = vout - vin;
rfft = fftn(vres);
ecut2 = get(mol,'ecut2');
gmask2 = FreqMask(mol,ecut2);
gm2  = get(gmask2,'gm');
gkk2 = get(gmask2,'gkk');
%
inz  = find(gm2~=0);
rfft(inz) = (0.8*gkk2(inz)./(gkk2(inz)+0.5)-0.8).*rfft(inz);
%
izr  = find(gkk2==0); 
rfft(izr) = -rfft(izr)/2.0;
%
vcor = ifftn(rfft);
vnew = vin + vcor + 0.8*vres;
