function e = getEntropy(occ,Tbeta);
%
% usage: e = getEntropy(occ, Tbeta);
%

%
% nonzero occupations
%
inz = find(occ > eps);
occnz = occ(inz);
e = sum(occnz.*log(occnz));
%
% zero occupations
%
iz = find(abs(occ) <= eps);
occz = occ(iz);
e = e + sum((1-occz).*log(1-occz));
e = e / Tbeta;
