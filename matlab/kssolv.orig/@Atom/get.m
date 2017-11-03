function val = get(a,attr_name)
%
% usage: val = get(a,attr_name);
%        retrieve various attributes of an Atom object a.
% e.g.  asymbol = get(a,'symbol') returns the atomic symbol associated with 
%       the Atom object a.
%       an = get(a,'anum') returns the atomic number associated with a
switch attr_name
  case 'symbol'
    val = a.symbol;
  case 'anum'
    val = a.anum;
  case 'amass'
    val = a.amass;
  case 'venum'
    val = a.venum;
  case 'iloc'
    val = a.iloc;
  case 'occs'
    val = a.occs;
  case 'occp'
    val = a.occp;
  case 'occd'
    val = a.occd;
  case 'iso'
    val = a.iso;
  case 'ic'
    val = a.ic;
  case 'isref'
    val = a.isref;
  case 'ipref'
    val = a.ipref;
  case 'idref'
    val = a.idref;
end;
