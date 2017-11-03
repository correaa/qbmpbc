function [vhart,vxc,uxc2,rho_out] = getvhxc(mol,rho_in)
%
% Usage: [vhart,vxc,uxc2,rho_out] = getvhxc(mol,rho_in)
%
% Purpose:
%    Compute the Hartree (vhart) and exchange correlation (vxc) 
%    potential energy based on the input charge density (rho_in)
%
% Input:
%  mol    --- a Molecule object
%  rho_in --- input charge density (3D)
%
% Output:
%  vhart   --- Hartree potential (3D)
%  vhxc    --- Exchange correlation potential (3D)
%  uxc2    --- Needed for Exchange correlation energy calculation (3D)
%  rho_out --- output charge density (3D, not that much different from input) 
%
% The Hartee potential is obtained by solving the Poisson equation
%                Laplacian*vhart = rho_in
% with periodic boundary condition in frequency space.
%
% The exchange correlation potential is computed by the Perdew & Zunger formula
%

%
[n1,n2,n3]=size(rho_in);
ecut2 = get(mol,'ecut2');
gmask2 = FreqMask(mol,ecut2);
gm2  = get(gmask2,'gm');
gkk2 = get(gmask2,'gkk');
%
[m1,m2,m3] = size(rho_in);
if (m1 ~= n1 | m2 ~= n2 | m3 ~= n3)
   error('get_hxc: dimension of the molecule does not match that of the rho');
end;
n123   = n1*n2*n3;
%
% Calculate Hartree potential
%
rhog   = fftn(rho_in);
w      = zeros(n1,n2,n3);
inz    = find(abs(gkk2) ~= 0);
w(inz) = 2*pi*rhog(inz)./gkk2(inz);
vhart  = real(ifftn(w));
%
% Calculate the exchange correlation potential
%
beta11  = 1.228383333;
beta22  = 0.44453333333;
b1   = -0.11673333;
c1   = 0.0026666667;
d1   = -0.0168;
gama = -0.2846;
beta1   = 1.0529;
beta2   = 0.3334;
tft  = 0.75;
cex  = -1.969490099;
a    = 0.0622;
b    = -0.096;
c    = 0.0040;
d    = -0.0232;

vc    = zeros(n1,n2,n3);
ec    = zeros(n1,n2,n3);
uxc2  = zeros(n1,n2,n3);
uxcca = zeros(n1,n2,n3);
%
rh3 = rho_in;
itny = find(rh3<1e-16);
if (~isempty(itny))
  rh3(itny) = 1e-16;
end; 
rho_out = rh3;
%
rh3   = rh3.^(1/3);
rs    = (pi*rho_out/tft).^(-1/3);
ibig  = find(rs >= 1);
isml  = find(rs <  1);

rsbig        = rs(ibig);
rootrs       = sqrt(rsbig);

vc(ibig)     = gama*(1+beta11*rootrs+beta22*rsbig)...
               ./(1+beta1*rootrs+beta2*rsbig).^2;

ec(ibig)     = gama./(1+beta1*rootrs+beta2*rsbig);

vexcor(ibig) = cex*rh3(ibig)+vc(ibig);

uxc2(ibig)   = (tft*cex*rh3(ibig)+ec(ibig))/2;
uxcca(ibig)  = vexcor(ibig)/2;

xlnrs        = log(rs(isml));
rssml        = rs(isml);

vc(isml)     = a*xlnrs+b1+c1*rssml.*xlnrs+d1*rssml;
ec(isml)     = a*xlnrs+b+c*rssml.*xlnrs+d*rssml;
vexcor(isml) = cex*rh3(isml)+vc(isml);
uxc2(isml)   = (tft*cex*rh3(isml)+ec(isml))/2;
uxcca(isml)  = vexcor(isml)/2;

vxc = uxcca;
