function KX = applyKIEP(H,X)
%
% Perform the multiplication of a Hamiltonian (H) and 
% wavefunction(s) psi.
%
% usage: KX = mtimes(H,X);
%
psi = get(X,'psi');
iscompact = get(X,'iscompact');
idxnz = get(X,'idxnz');
m1=get(X,'n1');
m2=get(X,'n2');
m3=get(X,'n3');
%
if (nargin == 2)
   gkk = H.gkk;
   gm  = H.gm;
   [n1,n2,n3] = size(H.vion);
   iz  = find(gm==0);
   if (iscell(psi))
      ncols = get(X,'ncols');
      for j = 1:ncols
         if (m1==n1 | m2==n2 | m3==n3)
            if (iscompact)
               psi3d = zeros(n1,n2,n3);
               psi3d(idxnz) = psi{j};
               psir  = ifftn(psi3d); 
            else
               psir  = ifftn(psi{j});
            end;
            vpsir = (H.vion+H.vext).*psir;
            if (iscompact)
               vpsi3d = fftn(vpsir);
               vpsi = vpsi3d(idxnz);
            else
               vpsi = fftn(vpsir);
               vpsi(iz) = 0;
            end;
            hpsi{j} = vpsi + gkk(idxnz).*psi{j};
            %
            % apply nonlocal pseudopotential
            %
            wqmat  = get(H,'wqmat');
            wqsign = get(H,'wqsign');
            if (~isempty(wqmat))
               hpsi{j} = hpsi{j} + wqmat*(wqsign.*(wqmat'*psi{j}));
            end;
         else
            error('The dimension of the Hamiltonian does not match that of the wave function');
         end;
      end;
   else
      [m1,m2,m3]=size(psi);
      if (m1==n1 | m2==n2 | m3==n3)
         psir  = ifftn(psi);
         vpsir = (H.vion+H.vext).*psir;
         vpsi  = fftn(vpsir);
         vpsi(iz) = 0;
         hpsi  = vpsi + gkk.*psi;
         %
         % apply nonlocal pseudopotential
         %
         wqmat  = get(H,'wqmat');
         wqsign = get(H,'wqsign');
         if (~isempty(wqmat))
            hpsi = hpsi + wqmat*(wqsign.*(wqmat'*psi));
         end;
      else
         error('The dimension of the Hamiltonian does not match that of the wave function');
      end;
   end;
   if (iscompact)
      KX = Wavefun(hpsi,n1,n2,n3,idxnz);   
   else
      KX = Wavefun(hpsi);   
   end;
else
   error('Ham: multiplication syntax: ham*wavefunction')
end;
