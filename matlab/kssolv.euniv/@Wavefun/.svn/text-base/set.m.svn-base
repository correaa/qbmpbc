function X = set(X,varargin)
%
% usage: X = set(X, attr_name1, value1, attr_name2, value2,...);
% set the attributes of a Wavefun object X.
%
Xargin = varargin;
while (length(Xargin)>=2),
   attr_name = Xargin{1};
   value     = Xargin{2};
   Xargin    = Xargin(3:end);
   switch attr_name
     case 'n1'
        X.n1 = value;
     case 'n2'
        X.n2 = value;
     case 'n3'
        X.n3 = value;
     case 'ncols'
        X.ncols = value;
     case 'nrows'
        X.ncols = value;
     case 'psi'
        X.psi = value;
     case 'idxnz'
        X.idxnz = value;
     case 'trans'
        X.trans = value;
     case 'iscompact'
        X.iscompact = value;
     otherwise
        fprintf('Attribute name: %s\n', attr_name);
        error('cannot be found');
   end;
end;
