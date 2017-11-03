function Y = subsref(X,S)
%
% usage: Y = X(S);
%
rows = S.subs{1};
cols = S.subs{2};
Ypsi = X.psi(rows,cols);
if (X.iscompact)
   % have to check this!!!!!, 
   % number of nonzeros may be reduced
   Y = Wavefun(Ypsi,X.n1,X.n2,X.n3,X.idxnz);
else
   Y = Wavefun(Ypsi);
end;
