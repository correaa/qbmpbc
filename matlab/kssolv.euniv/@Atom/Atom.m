function a = Atom(varargin)
%
% Atom constructor
%
switch (nargin)
  case 0
     a.symbol = [];
     a.anum   = [];
     a.amass  = [];
     a.venum  = [];
     a.iloc   = [];
     a.occs   = [];
     a.occp   = [];
     a.occd   = [];
     a.iso    = [];
     a.ic     = [];
     a.isref  = [];
     a.ipref  = [];
     a.idref  = [];
     a = class(a, 'Atom');
  case 1
     a.symbol = [];
     a.anum   = [];
     a.amass  = [];
     a.venum  = [];
     a.iloc   = [];
     a.occs   = [];
     a.occp   = [];
     a.occd   = [];
     a.iso    = [];
     a.ic     = [];
     a.isref  = [];
     a.ipref  = [];
     a.idref  = [];
     if ( ischar(varargin{1}) )
        a.symbol = varargin{1};
        [a.anum,a.amass] = alookup(varargin{1});
     else
        a.symbol = slookup(varargin{1});
        [a.anum,a.amass] = alookup(a.symbol);
     end; 
     [venum,iloc,occs,occp,occd,iso,ic,isref,ipref,idref] = elookup(varargin{1});
     a.venum = venum;
     a.iloc  = iloc; 
     a.occs  = occs;
     a.occp  = occp;
     a.occd  = occd;
     a.iso   = iso;
     a.ic    = ic;
     a.isref = isref;
     a.ipref = ipref;
     a.idref = idref;
     a = class(a,'Atom'); 
  otherwise
     error('Cannot have more than one arguement');
end;
