function rho = getcharge(mol,X,occ)
%
% Usage: rho = getcharge(mol,X,occ);
%
% Purpose:
%    Compute charge density from a given set of wavefunctions.
%
% Input:
%    mol --- a Molecule object
%     X  --- a set of wave functions (Wavefun object)
%  nspin --- spin type, if nspin=1, there is no distinction between up and down
%            spin. if npsin=2, up and down spins have to be treated differently.
%
% Ouput:
%    rho --- charge density (3D array)
%

if (nargin == 1)
   %
   %  construct rho from the atomic charges
   %
   pseudovar = pseudoinit(mol);
   [vion,vionT,rho] = getvion(mol,pseudovar);
   %
else
   vol = get(mol,'vol');
   psi = get(X,'psi');
   nspin = get(mol,'nspin');
   nocc  = getnel(mol)/2*nspin;
   n1 = get(X,'n1');
   n2 = get(X,'n2');
   n3 = get(X,'n3');
   n123 = n1*n2*n3;
   rho  = zeros(n1,n2,n3);
   iscompact = get(X,'iscompact');
   idxnz     = get(X,'idxnz');
   if (nargin == 2)
      %
      %  fully occupy the leading nocc states
      %
      for j = 1:nocc
         if (iscompact)
            psi3d = zeros(n1,n2,n3);
            psi3d(idxnz) = psi{j};
            psir = ifftn(psi3d)*n123;
         else
            psir = ifftn(psi{j})*n123;
         end;
         rho = rho + (2.0/nspin)*abs(psir).^2;
      end;
      rho = rho/vol;
   elseif (nargin == 3)
      %
      %  use partial occupation numbers stored in occ
      %
      ncols = get(X,'ncols');
      for j = 1:ncols
         if (iscompact)
            psi3d = zeros(n1,n2,n3);
            psi3d(idxnz) = psi{j};
            psir = ifftn(psi3d)*n123;
         else
            psir = ifftn(psi{j})*n123;
         end;
         rho = rho + occ(j)*(2.0/nspin)*abs(psir).^2;
      end;
      rho = rho/vol;
   else
      error('getcharge usage: getcharge(mol,X,occ);');
   end;
end;
