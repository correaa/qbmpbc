function [U,S,V]=svd(X,varargin)
%
% usage: [U,S,V] =svd(X);
%        [U,S,V] =svd(X,0);
%
psi = get(X,'psi');
n1  = get(X,'n1');
n2  = get(X,'n2');
n3  = get(X,'n3');
ncols = size(psi,2);
iscompact = get(X,'iscompact');
idxnz = get(X,'idxnz');

n123 = n1*n2*n3;
trans = get(X,'trans');

if (iscompact) 
   	Xmat = zeros(X.nrows,X.ncols);
     for j = 1:X.ncols
        Xmat(:,j) = X.psi{j};
     end;
  else
     Xmat = zeros(n123,X.ncols);
     for j = 1:X.ncols
        Xmat(:,j) = reshape(X.psi{j},n123,1);
     end;
  end
%
if (nargin == 1)
   [Umat,S,V] = svd(Xmat); 
elseif (nargin == 2)
   opt = varargin{1};
   if (opt==0) 
      [Umat,S,V]=svd(Xmat,0);
   else
      error('Invalid option for svd\n');
   end;
   if (iscompact)
      for j = 1:ncols
         psi{j} = Umat(:,j);
      end;
      U = Wavefun(psi,n1,n2,n3,idxnz);
   else
      for j = 1:ncols
         psi{j} = reshape(Umat(:,j), n1, n2, n3);
      end;
      U = Wavefun(psi);
   end;
else
   error('Incorrect number of arguements\n'); 
end;
