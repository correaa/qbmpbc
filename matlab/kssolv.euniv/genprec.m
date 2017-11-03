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
mpbc = get(mol,'mpbc');

if mpbc == 0
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

else
gmask = FreqMask(mol);
gm  = get(gmask,'gm');
%gkk = get(gmask,'gkk');
idxnz = find(gm~=0);
n1 = get(mol,'n1');
n2 = get(mol,'n2');
n3 = get(mol,'n3');
C  = get(mol,'supercell');
CI = inv(C);
%
gXhat = 2*pi*CI(1,1)*[0:n1*n2/2-1, -n1*n2/2:-1]/n2;
if n3 == 1
    gZ = 0;
else
    gZ = 2*pi*CI(3,3)*[0:n3/2-1,-n3/2:-1];
end
gkk = zeros(n1*n2*n3,1);
l=1;
for j=1:n3
    for i=1:n1*n2
        gkk(idxnz(l))=gXhat(i)^2+gZ(j)^2;
        l=l+1;
    end
end
Ek = 0.5;
X  = gkk/Ek;
Y  = 27.0 + X.*(18.0 + X.*(12.0 + 8.0*X));
p  = Y./(Y + 16.0*X.^4); 
prec = Wavefun({p},n1,n2,n3,idxnz);
end