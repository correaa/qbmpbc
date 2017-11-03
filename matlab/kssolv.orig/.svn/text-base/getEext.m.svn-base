function Eext = getEext(mol,rho,vext)
%
% Usage:
%    Eext = getEext(rho, vext);
%
% Purpose:
%    Compute the potential energy induced by an external potential vext.
%
% Input:
%  mol  --- Molecule object 
%  rho  --- input charge density (3D array)
%  vext --- external potential (3D array)
%
% Output:
%  Eext --- Potential energy induced by vext (scalar)
%

%
[n1,n2,n3]= size(rho);
vol       = get(mol,'vol');
Eext      = sum3d(rho.*vext)*vol/(n1*n2*n3);
Eext      = real(Eext);
