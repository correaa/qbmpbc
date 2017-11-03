function Y = ctranspose(X)
%
% Y = -X;
%
% usage: Y = -X;
%
if (nargin == 1)
   psi = X.psi;
   ncol = size(psi,2);
   for j = 1:ncol
     psi{j} = -psi{j};
   end 
   if (X.iscompact)
      idxnz = X.idxnz;
      n1 = X.n1;
      n2 = X.n2;
      n3 = X.n3;
      Y = Wavefun(psi,n1,n2,n3,idxnz);
   else
      Y = Wavefun(psi);
   end;
else
   error('Wavefun: complex transpose takes only one arguement')
end;
