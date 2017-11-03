function prec = genprec(mol)
%
% Usage: prec = genprec(mol);
%
% Purpose: 
%   Generate a preconditioner for LOBPCG and DCM in the frequency space
%   using the Teter recipe
%
% Input:
%   mol  --- a Molecule object
%   prec --- preconditioner stored as a Wavefun object.
%
gmask = FreqMask(mol);
gm  = get(gmask,'gm');
gkk = get(gmask,'gkk');
idxnz = find(gm~=0);
n1 = get(mol,'n1');
n2 = get(mol,'n2');
n3 = get(mol,'n3');
%
Ek = 0.5;
X  = gkk(idxnz)/Ek;
Y  = 27.0 + X.*(18.0 + X.*(12.0 + 8.0*X));
p  = Y./(Y + 16.0*X.^4); 
prec = Wavefun({p},n1,n2,n3,idxnz);
