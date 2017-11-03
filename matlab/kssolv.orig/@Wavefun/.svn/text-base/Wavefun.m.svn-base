function X = Wavefun(varargin)
%
% Wavefuniltonian constructor
%
% usage:  X = Wavefun();            constructs a null Wavefun object
%         X = Wavefun(n1,n2,n3);    constructs a single wave function
%         X = Wavefun(n1,n2,n3,ncols);  constructs ncols wave functions
%         X = Wavefun(psi);         psi is a cell array of 3-D arrays
%                                   or simply a 3-D array 
%         X = Wavefun(psi,n1,n2,n3,idxnz); 
%                                   constructs a Wavefun using a compact
%                                   storage scheme
%
switch (nargin)
  case 0
     X.n1    = [];  % number of sampling point in x,y,z
     X.n2    = [];
     X.n3    = []; 
     X.nrows = [];  % number of elements (rows) in each wavefunctions
                    % When a compact representation is used, nrows is the number 
                    % of nonzero frequency vectors; When a non-compact 
                    % representation 
                    % is used, nrows is undefined.
     X.ncols = [];  %number of wavefunctions 
     X.psi   = [];  %store the Fourier coefficients
     X.idxnz = [];  %location of the nonzeros in 3-D array if
                    %compac storage is used
     X.trans = 0;   % 0 - no transpose, 1 - conjugate transpose
     X.iscompact=0; % default to no compact storage
     X = class(X, 'Wavefun');
  case 1
     %
     % the input argument is 3-D Fourier coefficients of wavefunctions
     %
     X.n1    = 0;
     X.n2    = 0;
     X.n3    = 0;
     X.ncols = 0;
     X.nrows = 0;
     X.idxnz = [];
     X.trans = 0;
     Xin = varargin{1}; 
     X.iscompact=0; % default to no compact storage
     if ( iscell(Xin) )
        X.psi  = Xin;
        X.ncols = size(Xin,2);
        [n1,n2,n3]=size(Xin{1});
        X.n1 = n1;
        X.n2 = n2;
        X.n3 = n3;
        X.nrows = [];
        X = class(X,'Wavefun'); 
     elseif ( isnumeric(Xin) )
        X.psi  = {Xin};
        X.ncols = 1;
        [n1,n2,n3]=size(Xin);
        X.n1 = n1;
        X.n2 = n2;
        X.n3 = n3;
        X.nrows = [];
        X = class(X,'Wavefun'); 
     else 
        error('The input must be a cell array or a 3-D array');
     end; 
  case 3
     %
     % create one null wavefunction of size n1xn2xn3.
     %
     X.n1  = varargin{1};
     X.n2  = varargin{2};
     X.n3  = varargin{3};
     X.ncols  = 1;
     X.nrows  = [];
     X.psi = [];
     X.idxnz = []; 
     X.trans = 0;
     X.iscompact=0; 
     X = class(X,'Wavefun'); 
  case 4
     %
     % create ncols null wavefunctions of size n1xn2xn3
     %
     X.n1  = varargin{1};
     X.n2  = varargin{2};
     X.n3  = varargin{3};
     X.ncols  = varargin{4};
     X.nrows  = [];
     X.psi = [];
     X.idxnz = []; 
     X.trans = 0;
     X.iscompact=0; 
     X = class(X,'Wavefun'); 
  case 5 
     %
     % create ncols 3-D wavefunctions of size n1xn2xn3 using
     % compact storage
     %
     X.nrows = [];
     ncols = size(varargin{1},2);  
     if (isnumeric(varargin{1}))
        psimat = varargin{1};
        X.nrows = size(psimat,1);
        for j = 1:ncols
          X.psi{j} = psimat(:,j); 
        end;
     elseif (iscell(varargin{1}))
        X.psi = varargin{1};
        X.nrows = size(X.psi{1},1);
     else
        error('wrong input data type');
     end;
     X.n1    = varargin{2};
     X.n2    = varargin{3};
     X.n3    = varargin{4};
     X.idxnz = varargin{5};
     X.ncols = ncols;

     X.trans = 0;
     X.iscompact=1; 
     X = class(X,'Wavefun'); 
  otherwise
     error('wrong number of arguements');
end;
