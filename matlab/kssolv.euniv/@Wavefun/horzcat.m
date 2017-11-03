function Y = horzcat(varargin)
%
% usage: Y = [X1 X2 X3];
% horizontal concatenation
%
m = length(varargin);
Ypsi = [];
for j = 1:m
   if ( isa(varargin{j},'Wavefun') )
      Ypsi = [Ypsi get(varargin{j},'psi')];
   else
      fprintf('All entries must be Wavefun objects\n');
   end;
end;
X1 = varargin{1};
if (~isempty(X1))
   if (X1.iscompact)
      Y = Wavefun(Ypsi,X1.n1,X1.n2,X1.n3,X1.idxnz); 
   else
      Y = Wavefun(Ypsi); 
   end;
end;
