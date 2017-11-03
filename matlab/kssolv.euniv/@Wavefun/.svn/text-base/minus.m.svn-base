function Z = minus(X,Y)
%
% Perform pointwise multiplication of two Wavefun objects
%              Z = X-Y;
%
if (nargin == 2)
   %
   % --- process the first operand ---
   %
   if ( isa(X,'Wavefun') & isa(Y,'Wavefun') )
      %
      % turn the cell arrays into a matrices
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
      end; 
      %
      Yn123 = Y.n1*Y.n2*Y.n3;
     	if (Y.iscompact)
          Ymat = zeros(Y.nrows,Y.ncols);
          for j = 1:Y.ncols
             Ymat(:,j) = Y.psi{j};
          end;
       else 
          Ymat = zeros(Yn123,Y.ncols);
          for j = 1:Y.ncols
             Ymat(:,j) = reshape(Y.psi{j},Yn123,1);
          end;
       end;
      if (X.trans==0 & Y.trans==0)
         Zmat = Xmat-Ymat;
      elseif (X.trans == 1 & Y.trans == 0)
         Zmat = Xmat'-Ymat;
      elseif (X.trans == 0 & Y.trans == 1)
         Zmat = Xmat-Ymat';
      else
         Zmat = Xmat'-Ymat';
      end;
      if (X.iscompact)  
         for j = 1:X.ncols
            Zpsi{j} = Zmat(:,j);
         end
         Z = Wavefun(Zpsi,X.n1,X.n2,X.n3,X.idxnz);
      else
         for j = 1:X.ncols
            Zpsi{j} = reshape(Zmat(:,j),X.n1,X.n2,X.n3);
         end
         Z = Wavefun(Zpsi);
      end;
   else
     error('The operands must be Wavefun objects');
   end;
else
   error('Wavefun: multiplication syntax: Z = X+Y');
end;
