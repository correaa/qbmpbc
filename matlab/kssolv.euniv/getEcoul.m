function Ecoul = getEcoul(mol,rho,vhart)
%
% Usage: Ecoul = getEcoul(mol, rho, vhart);
%
% Purpose:
%    Compute the Coulomb potential energy induced by vhart
%
% Input:
%  mol    --- a Molecule object
%  rho    --- input charge density (3D array)
%  vhart  --- Hartee (Coulomb) potential (3D array)
%
% Output:
%  Ecoul  --- Coulomb potential energy induced by vhart.
%

[n1,n2,n3] = size(rho);
vol        = get(mol,'vol');
Ecoul      = sum3d(rho.*vhart)/2;
Ecoul      = Ecoul*vol/(n1*n2*n3);
Ecoul      = real(Ecoul);
