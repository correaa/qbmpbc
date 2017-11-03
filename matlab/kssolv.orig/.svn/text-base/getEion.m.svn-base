function Eion = getEion(mol,rho,vion)
%
% Usage: Eion = getEion(rho, vion);
%
% Purpose:
%    Compute the ionic potential energy
%
% Input:
%    mol  --- a Molecule object
%    rho  --- input charge density (3D array)
%    vion --- ionic potential (3D array)
%
% Output:
%   Eion  --- Ionic potential energy (scalar)
%

[n1,n2,n3]= size(rho);
vol       = get(mol,'vol');
Eion      = sum3d(rho.*vion)*vol/(n1*n2*n3);
Eion      = real(Eion);
