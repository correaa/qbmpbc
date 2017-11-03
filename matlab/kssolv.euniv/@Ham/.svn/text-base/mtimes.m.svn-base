function HX = mtimes(H,X)
%
% Perform the multiplication of a Hamiltonian (H) and 
% wavefunction(s) X.
%
% usage: HX = mtimes(H,X);
%
if (nargin == 2 && H.mpbc == 0)
   disp('zero magnetic field');
   gkk = H.gkk;
   gm  = H.gm;
   idxnz = find(gm~=0);
   [n1,n2,n3] = size(H.vtot);
   wqmat = [];
   wqsign = [];
   if ( isa(X,'Wavefun') & get(X,'iscompact') )
      ncols = get(X,'ncols');
      Xpsi = get(X,'psi');
      for j = 1:ncols
         m1 = get(X,'n1');
         m2 = get(X,'n2');
         m3 = get(X,'n3');
         if (m1==n1 & m2==n2 & m3==n3)
            Xpsi3d = zeros(n1,n2,n3);
            Xpsi3d(idxnz) = Xpsi{j};
            Xr3d   = ifftn(Xpsi3d);
            vXr3d  = H.vtot.*Xr3d;
            vX3d   = fftn(vXr3d);
            vX     = vX3d(idxnz);
            HXpsi{j} = vX + gkk(idxnz).*Xpsi{j};
            %
            % apply nonlocal pseudopotential
            %
            wqmat  = get(H,'wqmat');
            wqsign = get(H,'wqsign');
            if (~isempty(wqmat))
               HXpsi{j} = HXpsi{j} + wqmat*(wqsign.*(wqmat'*Xpsi{j}));
            end;
         else
            error('The dimension of the Hamiltonian does not match that of the wave function');
         end;
      end;
      HX = Wavefun(HXpsi,n1,n2,n3,idxnz);
   else
      error('The Hamiltonian must operate on a Wavefun object');
   end;
elseif (nargin == 2 && H.mpbc > 0)
    disp('nonzero magnetic field');
    HX = MPBC_mtimes(H,X);
else
   error('Ham: multiplication syntax: Y=H*X')
end;
