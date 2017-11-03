function VX = applyNP(H,X)
%
% Perform the multiplication of a Hamiltonian (H) and 
% wavefunction(s) X.
%
% usage: VX = applyNP(H,X);
%
if (nargin == 2)
   gkk = H.gkk;
   gm  = H.gm;
   [n1,n2,n3] = size(H.vnp);
   iz  = find(gm==0);
   if (isa(X,'Wavefun'))
      ncols = get(X,'ncols');
      Xpsi  = get(X,'psi');
      iscompact = get(X,'iscompact');
      idxnz     = get(X,'idxnz');
      %
      for j = 1:ncols
         m1 = get(X,'n1');
         m2 = get(X,'n2');
         m3 = get(X,'n3');
         if (m1==n1 & m2==n2 & m3==n3)
            if (iscompact)
               Xpsi3d = zeros(n1,n2,n3);
               Xpsi3d(idxnz) = Xpsi{j};
               Xr = ifftn(Xpsi3d); 
            else 
               Xr = ifftn(Xpsi{j});
            end;
            VXr  = (H.vnp).*Xr;
            if (iscompact)
               VXpsi3d = fftn(VXr);
               VXpsi{j} = VXpsi3d(idxnz);
            else
               VXpsi{j} = fftn(VXr);
               VXpsi{j}(iz) = 0;
            end;
         else
            error('The dimension of the Hamiltonian does not match that of the wave function');
         end;
      end;
      if (iscompact)
         VX = Wavefun(VXpsi,n1,n2,n3,idxnz);
      else
         VX = Wavefun(VXpsi);
      end;
   else
      error('The Hamiltonian must operate on a Wavefun object');
   end;
else
   error('Ham: applyNP syntax: VX=applyNP(H,X)');
end;
