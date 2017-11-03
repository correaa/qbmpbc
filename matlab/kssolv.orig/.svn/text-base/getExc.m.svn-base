function Exc = getExc(mol,rho,uxc2)
%
% Usage: Exc = getExc(mol,rho,vexc);
%
% Purpose:
%    Compute the exchange correlation energy
%
% Input:
%  mol  --- Molecule object
%  rho  --- input charge density (3D array)
%  vexc --- exchange correlation potential (3D array)
%
% Output:
%  Exc  --- Exchange correlation energy (scalar)
%

[n1,n2,n3]= size(rho);
n123      = n1*n2*n3;
vol       = get(mol,'vol');
Exc       = sum3d(uxc2.*rho)*vol/n123;
