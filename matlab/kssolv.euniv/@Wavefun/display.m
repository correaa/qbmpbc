function display(X)
%
% display method for atom
%
fprintf('Wavefun object\n');
fprintf('number of wave functions: ncols = %d\n', X.ncols);  
fprintf('wave function size: n1 = %d, n2 = %d, n3 = %d\n', X.n1, X.n2, X.n3);
if (X.trans)
   fprintf('transpose: Yes\n');
else
   fprintf('transpose: No\n');
end;
%
if (X.iscompact)
   fprintf('compact storage: Yes\n');
   fprintf('nrows = %d, ncols = %d\n', X.nrows, X.ncols);
else
   fprintf('compact storage: No\n');
end;
