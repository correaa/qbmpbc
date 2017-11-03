function H = set(H,varargin)
%
% usage: H = set(H, attr_name1, value1, attr_name2, value2,...);
% set the attributes of a Hamiltonian object H.
%
Hargin = varargin;
while (length(Hargin)>=2),
   attr_name = Hargin{1};
   value     = Hargin{2};
   Hargin    = Hargin(3:end);
   switch attr_name
     case 'gmask'
        H.gkk = get(value,'gkk');
        H.gm  = get(value,'gm'); 
     case 'vion'
        H.vion = value; 
     case 'vext'           
        H.vext = value;    
     case 'vtot'
        H.vtot = value;
     case 'vnp'
        H.vnp = value;
     case 'rho'
        H.rho = value;
        % 
        % Calculate Hartree and exchange-correlation potential
        % from the known density rho
        % 
        [vhart,vxc,uxc2,rho]=getvhxc(mol,abs(H.rho));
        %
        % Save a copy of the Hartree and exchange-correlation
        % separately for DCM updating purposes
        %
        H.vnp = vhart+vxc;
        %
        % vtot does not include non-local ionici potential
        % 
        H.vtot = getvtot(mol, H.vion, H.vext, vhart, vxc);
    case 'mpbc'
       H.mpbc = value;
    case 'mpbc_denom'
       H.mpbc_denom = value;
    case 'supercell'
       H.supercell = value;
    otherwise
        fprintf('Attribute name: %s\n', attr_name);
        error('cannot be found');
   end;
end;
