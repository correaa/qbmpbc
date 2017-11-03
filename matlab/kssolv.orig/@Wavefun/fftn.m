function Y = fftn(X)
%
% Usage: Y = fftn(X)
%
% Purpose:
%    Perform multidimensional FFT on a Wavefun object. 
%
% Input:
%    X  --- a Wavefun object, may or may not be stored using
%           a compact format that keeps only the Fourier coefficients
%           associated with frequency vectors that satisfy the 
%           kinetic energy cutoff.
% Output:
%
%    Y  --- a Wavefun object. Y is not represented by a compact storage format.
%
xpsi = get(X,'psi');
ncols = get(X,'ncols');
fpsi = [];
iscompact = get(X,'iscompact');
if (iscompact)
   idxnz = get(X,'idxnz');
   n1 = get(X,'n1');
   n2 = get(X,'n2');
   n3 = get(X,'n3');
   for j = 1:ncols
     psi3d = zeros(n1,n2,n3);
     psi3d(idxnz) = xpsi{j};
     fpsi{j} = fftn(psi3d);
   end; 
else
   for j = 1:ncols
     fpsi{j} = fftn(xpsi{j});
   end;
end;
%Y = fpsi;
Y = Wavefun(fpsi);
