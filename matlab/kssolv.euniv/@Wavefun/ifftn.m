function Y = ifftn(X)
%
% usage: beta = ifftn(X)
%
psi  = get(X,'psi');
ncols  = get(X,'ncols');
iscompact = get(X,'iscompact');
if (iscompact) 
   n1 = get(X,'n1');
   n2 = get(X,'n2');
   n3 = get(X,'n3');
   idxnz = get(X,'idxnz');
   for j = 1:ncols
     psi3d = zeros(n1,n2,n3);
     psi3d(idxnz) = psi{j};
     fpsi{j} = ifftn(psi3d);
   end;
else
   for j = 1:ncols
     fpsi{j} = ifftn(psi{j});
   end;
end;
Y = Wavefun(fpsi);
