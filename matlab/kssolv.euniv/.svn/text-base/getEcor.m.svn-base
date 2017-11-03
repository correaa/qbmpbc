function Ecor = getEcor(mol, rho, vtot, vion, vext)
%
% Usage: Ecor = getEcor(mol, rho, vtot, vion, vext)
%
% Purpose:
%    Computes an energy correction term
%
% Input:
%    mol  --- a Molecule object
%    rho  --- charge density (3D array)
%    vtot --- total local potential (3D array)
%    vion --- local ionic potential (3D array)
%    vext --- external potential (3D array)
%
% Ouptut:
%    Ecor --- correction energy
%

[n1,n2,n3]= size(rho);
vol       = get(mol,'vol');
Ecor      = sum3d( (vion+vext-vtot).*rho );
Ecor      = Ecor*vol/(n1*n2*n3);
Ecor      = real(Ecor);
