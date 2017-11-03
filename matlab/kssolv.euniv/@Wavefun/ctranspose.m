function Y = ctranspose(X)
%
% Complex conjugate transpose
%
% usage: Y = X';
%
if (nargin == 1)
   psi = X.psi;
   if (X.iscompact)
      idxnz = X.idxnz;
      n1 = X.n1;
      n2 = X.n2;
      n3 = X.n3;
      Y = Wavefun(psi,n1,n2,n3,idxnz);
   else
      Y = Wavefun(psi);
   end;
   Y.trans = 1;
else
   error('Wavefun: complex transpose takes only one arguement')
end;
