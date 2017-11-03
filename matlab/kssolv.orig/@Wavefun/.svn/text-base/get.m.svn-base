function val = get(X,attr_name)
%
% usage: val = get(X,attr_name);
%        retrieve various attributes of a Wavefun object X.
% Wave functions are stored as a cell array of 3-D arrays.
% e.g.  n1  = get(X,'n1') returns the first dimension of the wavefunction
%              in the Xamiltonian.
%       psi = get(X,'psi') returns the actual cell array that
%             contains 3-D wave functions.
switch attr_name
  case 'n1'
    val = X.n1;
  case 'n2'
    val = X.n2;
  case 'n3'
    val = X.n3;
  case 'ncols'
    val = X.ncols;
  case 'nrows'
    val = X.nrows;
  case 'psi'
    val = X.psi;
  case 'idxnz'
    val = X.idxnz;
  case 'trans'
    val = X.trans;
  case 'iscompact'
    val = X.iscompact;
  otherwise
    error('invalid attribute');
end;
