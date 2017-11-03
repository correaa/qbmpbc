function Xabs = abs(X)
%
% usage: Xabs = abs(X)
%
for j = 1:X.k
   apsi{j} = abs(X.psi{j});
end;
if (X.iscompact)
   Xabs = Wavefun(apsi,X.n1,X.n2,X.n3,X.idxnz);
else
   Xabs = Wavefun(apsi);
end;
