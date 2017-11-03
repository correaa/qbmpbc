function vtot = getvtot(mol,vion,vext,vhart,vxc);
%
% Usage: vtot = getvtot(mol,vion,vext,vhart,vxc);
%
% Purpose:
%   Compute the total potential energy and keep it band limited
%
% Input: 
%     mol   a Molecular object
%    vion   a 3D array that contains the local ionic potential
%    vext   a 3D array that contains external potential. In most
%           case, this should be an empty array
%    vhart  a 3D array that contains the Hartree (Coulomb) potential.
%    vxc    a 3D array that contains the exchange-correlation potential
%
% Output:
%    vtot   a 3D array that contains the total (local) potential
%

[n1,n2,n3]=size(vion);
ecut2  = get(mol,'ecut2');
gmask2 = FreqMask(mol,ecut2);
gm2  = get(gmask2,'gm');
vtot = vion + vhart + vxc;
if (~isempty(vext))
  vtot = vtot + vext;
end;
vsum = sum3d(vtot);
vsum = vsum/(n1*n2*n3);  

% keep vtot band limited (sinc interpolation)
iz = find(gm2==0);
vfft = fftn(vtot);
vfft(iz) = 0;
vtot = ifftn(vfft);
vtot = vtot - vsum;
