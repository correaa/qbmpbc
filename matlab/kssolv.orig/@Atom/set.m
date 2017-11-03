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
     case 'symbol'
        if ( ischar(value) )
           a.symbol = value;
           a.anum   = alookup(value);
           [venum,iloc,occs,occp,occd,iso,ic,isref,ipref,idref]=elookup(value);
           a.venum = venum;
        else
           error('the symbol attribute must be followed by a char string arguement');
        end;
     case 'anum'
        if ( isnumeric(value) )
           a.anum   = value;
           a.symbol = slookup(value);
           [venum,iloc,occs,occp,occd,iso,ic,isref,ipref,idref]=elookup(value);
           a.venum  = venum;
        else
           error('the anum attribute must be followed by a numeric argument');
        end;
% we should not allow venum to be set by a user 
%     case 'venum'
%        if ( isnumeric(value) )
%           a.venum = value;
%        else
%           error('the venum attribute must be followed by a numeric argument');
%        end;
     otherwise
        fprintf('Attribute name: %s\n', attr_name);
        error('cannot be found');
   end;
end;
