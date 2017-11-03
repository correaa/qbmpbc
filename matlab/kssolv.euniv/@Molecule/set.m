function a = set(a,varargin)
%
% usage: a = set(a, attr_name1, value1, attr_name2, value2,...);
% set the attributes of an Atom object a.
% When the atomic symbol is set by the user, the method will set
% the corresponding atomic number automatically and vice versa
%
atom_argin = varargin;
while (length(atom_argin)>=2),
   attr_name = atom_argin{1};
   value     = atom_argin{2};
   atom_argin = atom_argin(3:end);
   switch attr_name
     case {'name', 'Name'}
        if ( ischar(value) )
           a.name = value;
        else
           error('the Name attribute must be followed by a char string arguement');
        end;
     case {'supercell'}
        if ( isnumeric(value) )
           [nrow,ncol]=size(value);
           if (nrow == 3 & ncol ==3)
              a.supercell  = value;
              a.vol       = abs(det(value));
           else
   	      error('supercell must be a 3 by 3 array');
           end;
        else
	   error('supercell must be a 3 by 3 array');
        end;
     case 'atomlist'
        a.atomlist = value;
        a.natoms   = length(value);
        a.nel      = getnel(a); 
        a.atypes   = getatypes(value);
     case 'xyzlist'
        %
        % check value to make sure it is a natoms by 3 array
        %
        a.xyzlist = value;
        a.natoms = size(value,1);
     case {'ecut','Ecut'}
        a.ecut  = value;
        a.ecut2 = 4*value;
        % we can estimate the minimum n1,n2,n3 values (based on
        % sampling theorem)
        %
        if ( ~isempty(a.supercell) )
           C = a.supercell; 
           fackpt = 2*sqrt(value)/pi;
           d1 = fackpt*norm(C(:,1));
           d2 = fackpt*norm(C(:,2));
           d3 = fackpt*norm(C(:,3));
           disp('resetting n1,n2,n3...');
           a.n1 = ceil(d1*0.5)*2;
           a.n2 = ceil(d2*0.5)*2;
           a.n3 = ceil(d3*0.5)*2;
           a.vext = zeros(a.n1,a.n2,a.n3);
        end;
     case {'ecut2', 'Ecut2'}
        a.ecut2 = value;
     case 'n1'
        % check against ecut
        a.n1 = value;
     case 'n2'
        % check against ecut
        a.n2 = value;
     case 'n3'
        % check against ecut
        a.n3 = value;
     case 'rcut'
        a.rcut = value;
     case {'nel','Nel'}
        a.nel = value;
     case {'vext','Vext'}
        a.vext = value;
     case {'nspin'}
        a.nspin = value;
     case {'temperature'}
        a.temperature = value;
     case 'mpbc'
        a.mpbc = value;
     case 'mpbc_denom'
        a.mpbc_denom = value;
     otherwise
        fprintf('Attribute name: %s\n', attr_name);
        error('cannot be found');
   end;
end;
