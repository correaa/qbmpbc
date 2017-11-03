function [Q,R]=qr(X,varargin)
%
% usage: Q    =qr(X,0);
%        [Q,R]=qr(X,0);
%        [Q,R]=qr(X);
%        Q    =qr(X);
%
psi = get(X,'psi');
n1  = get(X,'n1');
n2  = get(X,'n2');
n3  = get(X,'n3');
ncols   = get(X,'ncols');
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
if (trans==1)
   Xmat = Xmat';
end;
%
if (nargin == 1)
   [Qmat,R] = qr(Xmat); 
elseif (nargin == 2)
   opt = varargin{1};
   if (opt==0) 
      [Qmat,R]=qr(Xmat,0);
   else
      error('Invalid option for qr\n');
   end;
   if (iscompact)
      for j = 1:ncols
         psi{j} = Qmat(:,j);
      end;
      Q = Wavefun(psi,n1,n2,n3,idxnz);
   else
      for j = 1:ncols
         psi{j} = reshape(Qmat(:,j), n1, n2, n3);
      end;
      Q = Wavefun(psi);
   end;
else
   error('Incorrect number of arguements\n'); 
end;
