function beta = norm(X,varargin)
%
% usage: beta = norm(X)
%        beta = norm(X,'fro');
%        beta = norm(X,1);
%
n123 = X.n1*X.n2*X.n3;
if (X.ncols > 0)
   beta = 0.0;
   if (X.iscompact)
      for j = 1:X.ncols
         beta = beta + norm(X.psi{j})^2;
      end;
   else
      for j = 1:X.ncols
         beta = beta + norm(reshape(X.psi{j},n123,1))^2;
      end;
   end;
   beta = sqrt(beta);
else
   error('input does not contain any wavefunction!');
end;
