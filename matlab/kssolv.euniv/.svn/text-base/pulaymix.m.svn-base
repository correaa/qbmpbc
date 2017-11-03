function [vnew,dvmat,vmat] = pulaymix(vin, vout, dvmat, vmat, iter, mixdim)
%
% usage: [vnew,dvmat,vmat] = pulaymix(vin, vout, dvmat, vmat, iter, mixdim);
%
% purpose: Pulay mixing for accelerating the SCF iteration
%
% arguements:
%  input: 
%     vin (3-D array) --- input potential
%    vout (3-D array) --- output potential produced from eigenvalue 
%                         calculation
%    dvmat (matrix)   --- work array that stores the previous vout - vin.
%                         Its dimension is the same as that of vmat.
%    vmat (matrix)    --- work array that stores the previous input potentials
%                         The dimension of vmat is n1*n2*n3 by mixdim,
%                         where [n1,n2,n3]=size(vin).
%    iter (integer)   --- current SCF iteration number
%    mixdim(integer)  --- the number of columns in dvmat and vmat. This is the 
%                         the mixing memory. Only information from the
%                         previous mixdim iterations are used.
%  output:
%    vnew (3-D array) --- new potential to be used for the next SCF iteration
%    dvmat   (matrix) --- updated work array that keeps the previous
%                         vout-vin
%    vmat   (matrix)  --- updated work array that keeps the previous 
%                         input potentials
%  

[n1,n2,n3]=size(vin);
n123 = n1*n2*n3;
%
% construct and modify vmat and dvmat
%
if (iter <= mixdim)
   vmat(:,iter) = reshape(vin,n123,1);
   dvmat(:,iter) = reshape(vout - vin,n123,1);
else
   % delete the first leading column
   vmat(:,1) = []; dvmat(:,1) = [];
   % append vin and vout-vin at the end of vmat and dvmat;
   vmat(:,mixdim) = reshape(vin,n123,1);
   dvmat(:,mixdim) = reshape(vout-vin,n123,1);   
end;

if (iter > 1)
   ibeg = 1;
   iend = min(iter,mixdim);
   % 
   B = dvmat(:,ibeg:iend)'*dvmat(:,ibeg:iend);
   condB = cond(B);
   %debug
   fprintf('cond(B) = %11.3e\n', condB);
   while (condB > 1e10 & ibeg < iend)
      ibeg = ibeg + 1; 
      B = dvmat(:,ibeg:iend)'*dvmat(:,ibeg:iend);
      condB = cond(B);
      % debug 
      fprintf('cond(B) = %11.3e\n', condB);
   end;
   ncols = iend-ibeg+1;
   A = [B ones(ncols,1); ones(1,ncols) 0];
   b = [zeros(ncols,1); 1];
   g = A\b;
   vnew = reshape(vmat(:,ibeg:iend)*g(1:ncols,1),n1,n2,n3); 
else
   vnew = vout;
end;
