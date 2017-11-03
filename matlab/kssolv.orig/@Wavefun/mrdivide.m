function Z = mrdivide(X,R)
%
% Perform Z = X/R, where is a matrix
% 
if (nargin == 2)
   %
   % --- process the first operand ---
   %
   if ( isa(X,'Wavefun') )
     %
      % turn the cell array into a matrix
      %
      Xn123 = X.n1*X.n2*X.n3;
      if (X.iscompact)
         Xmat = zeros(X.nrows,X.ncols);
         for j = 1:X.ncols
            Xmat(:,j) = X.psi{j};
         end;
      else
         Xmat = zeros(Xn123,X.ncols);
         for j = 1:X.ncols
            Xmat(:,j) = reshape(X.psi{j},Xn123,1);
         end;
      end
      if (isnumeric(R))
         if (X.trans==0)
            Zmat = Xmat/R;
         else
            Zmat = Xmat'/R;
         end;
         %
         % --- convert to Wavefun object
         % 
         nz = size(Zmat,2);
         if (X.iscompact)
            for j = 1:nz
               Zpsi{j} = Zmat(:,j);
            end;
            Z = Wavefun(Zpsi,X.n1,X.n2,X.n3,X.idxnz); 
         else  
            for j = 1:nz
               Zpsi{j} = reshape(Zmat(:,j), X.n1, X.n2, X.n3);
            end;
            Z = Wavefun(Zpsi); 
         end; 
      else
         error('The second operand must be a matrix');
      end;
   else
     error('The first operand must be a Wavefun object');
   end;
else
   error('Wavefun: multiplication syntax: Z = X/R');
end;
