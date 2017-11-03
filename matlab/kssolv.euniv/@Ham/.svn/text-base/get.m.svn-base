function val = get(H,attr_name)
%
% usage: val = get(H,attr_name);
%        retrieve various attributes of an Hamiltonian object H.
% e.g.  vion = get(H,'vion') returns the ionic potential stored
%              in the Hamiltonian.
%       vtot = get(H,'vtot') returns the aggregated potential stored in 
%              the Hamiltonian 
switch attr_name
  case 'gkk'
    val = H.gkk;
  case 'gm'
    val = H.gm;
  case 'gkk2'
    val = H.gkk2;
  case 'gm2'
    val = H.gm2;
  case 'vion'
    val = H.vion;
  case 'wqmat'
    val = H.wqmat;
  case 'wqsign'
    val = H.wqsign;
  case 'vext'      % bhlee
    val = H.vext;  % bhlee
  case 'vtot'
    val = H.vtot;
  case 'rho'
    val = H.rho;
  case 'Ehxc'
    val = H.Ehxc;
  case 'mpbc'
    val = H.mpbc;
  case 'mpbc_denom'
    val = H.mpbc_denom;
  case 'supercell'
    val = H.supercell;
  otherwise
    error('invalid attribute');
end;
