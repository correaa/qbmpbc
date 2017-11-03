function X = subsasgn(X,S,Y)
%
% usage: X(index) = Y;
%
rows = S.subs{1};
cols = S.subs{2};
X.psi(rows,cols) = Y.psi;
if (X.iscompact & Y.iscompact)
   X = Wavefun(X.psi,X.n1,X.n2,X.n3,X.idxnz);
else
   X = Wavefun(X.psi);
end;
