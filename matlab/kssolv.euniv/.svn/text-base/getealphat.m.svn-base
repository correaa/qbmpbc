function Ealphat = getealphat(mol)
%
% Calculate the Ealpha energy
%
% usage: Ealphat = getealphat(mol);
%        where mol is a Molecule object;
%

%
% check the validity of the input argument.
%
%...

%fprintf('Calculating Ealphat energy...\n');
%t0      = cputime;
Ealphat = 0;
%
% the number of different atom types
%
elements = zeros(110,1);
ealpha   = zeros(110,1);
%
atomlist = get(mol,'atomlist');
natoms   = get(mol,'natoms');
vol      = get(mol,'vol');
% 
for j = 1:natoms
   a = atomlist(j);
   anum = get(a,'anum');
   elements(anum) = 1;
end;
ntypes = sum(elements);

inz = find(elements==1);
for j = 1:ntypes;
   anum = inz(j);
   a = Atom(anum);
   iloc = get(a,'iloc');
   vwr  = loadppt(a); 
   nrr = size(vwr,1);
   r   = vwr(:,1);
   vs  = vwr(:,2);
   vp  = vwr(:,3);
   vd  = vwr(:,4);
   ws  = vwr(:,5);
   wp  = vwr(:,6);
   wd  = vwr(:,7);
   %
   s = 0.0;
   ch = get(a,'venum');
   if (iloc == 1) 
      vloc = vs; 
   elseif (iloc == 2)
      vloc = vp;
   end;
   %for i = 2:nrr-1
   %   if (r(i) < 15.0)
   %      s=s+(ch*r(i)+vloc(i)*r(i)^2)*(r(i+1)-r(i-1))/2;
   %   end;
   %end;
   %
   % vectorize the above loop below
   %
   i15 = find(r(2:nrr-1)<15.0)+1;
   s = sum( (ch*r(i15)+vloc(i15).*r(i15).^2).*(r(i15+1)-r(i15-1))/2 ); 
   %
   ealpha(inz(j)) = s*4*pi;
   clear a;
end;
%
ch = 0.0;
for j = 1:natoms
   a = atomlist(j);
   anum = get(a,'anum');
   Ealphat = Ealphat + ealpha(anum);
   venum = get(a,'venum');
   ch = ch + venum;
end;
Ealphat = Ealphat*ch/vol;
%fprintf('Time = %11.3e\n', cputime - t0);
