function [wqmat,wqsign] = getwq(mol,pseudovar)
%
% Usage: [wqmat,wqsign] = getwq(mol,pseudovar);
%
% Purpose:
%   Computes the Kleiman-Bylander nonlocal ionic potential 
%   from atomic pseudopotential
%
% Input:
%    mol       --- an Molecule object
%    pseudovar --- pseudopotential information
%                  often obtained from pseudovar = pseudoinit(mol);
% Output:
%    wqmat     --- an ng x (9*natoms) array 
%    wqsign    --- an 9*natoms vector of 1's and -1's
%
fprintf('Calculating non-local ionic potential...\n');
t0    = cputime;
wqmat  = [];
wqsign = [];
% 
n1       = get(mol,'n1');
n2       = get(mol,'n2');
n3       = get(mol,'n3');
vol      = get(mol,'vol');
atomlist = get(mol,'atomlist');
natoms   = get(mol,'natoms');
vol      = get(mol,'vol'); 
%
mnq    = 4000;

if (~isempty(atomlist))
   xyzlist  = get(mol,'xyzlist');
   %
   gmask = FreqMask(mol);
   ng    = get(gmask,'ng');
   gm    = get(gmask,'gm');
   gkk   = get(gmask,'gkk');
   gkx   = get(gmask,'gkx');
   gky   = get(gmask,'gky');
   gkz   = get(gmask,'gkz');
   inz   = find(gm~=0);
   %
   wq2  = pseudovar.wq2;
   qi3  = pseudovar.qi3';
   isNL = pseudovar.isNL; 
   %
   wqmask = zeros(ng,9,natoms); 
   %
   % figure out the number of different types of atoms
   % sort the atomic numbers in ascending order
   %
   elements = zeros(110,1);
   for j = 1:natoms   
      a = atomlist(j);  
      anum = get(a,'anum'); 
      elements(anum) = 1; 
   end; 
   atypes = find(elements==1);
   %
   % atypes should list all types ordered by atomic numbers
   %
   for ia = 1:natoms
      a = atomlist(ia);
      x = xyzlist(ia,1);
      y = xyzlist(ia,2);
      z = xyzlist(ia,3);
      anum = get(a,'anum');
      %
      % identify the type the atom
      % 
      itype = find(atypes==anum);
      %fprintf('ia = %d, anum = %d, itype = %d\n', ia, anum, itype);
      % 
      isref = get(a,'isref');
      ipref = get(a,'ipref');
      idref = get(a,'idref');
      %
      iloc = get(a,'iloc');
      if (iloc == 1)
         isref = 0;
      elseif (iloc ==2)
         ipref = 0;
      else
         idref = 0;
      end; 
      %
      % phase vector
      %
      phvec = exp((gkx(inz)*x + gky(inz)*y + gkz(inz)*z)*sqrt(-1)); 
      %
      qvec = sqrt(gkx(inz).^2+gky(inz).^2+gkz(inz).^2);
      %
      iqvec = fix(1 + qvec*(mnq-1)/qi3(mnq));
      %
      xvec = (qvec-qi3(iqvec))./(qi3(iqvec+1)-qi3(iqvec));
      %
      f1vec = 1-xvec-0.5*xvec.*(1-xvec);
      f2vec = xvec+xvec.*(1-xvec);
      f3vec = -0.5*xvec.*(1-xvec);
      %
      ysvec = wq2(iqvec,1,itype).*f1vec + wq2(iqvec+1,1,itype).*f2vec ...
            + wq2(iqvec+2,1,itype).*f3vec;
      ypvec = wq2(iqvec,2,itype).*f1vec + wq2(iqvec+1,2,itype).*f2vec ...
            + wq2(iqvec+2,2,itype).*f3vec;
      ydvec = wq2(iqvec,3,itype).*f1vec + wq2(iqvec+1,3,itype).*f2vec ...
            + wq2(iqvec+2,3,itype).*f3vec;

      iqs = find(qvec<1e-6);
      if (isref == 1)
         wqmask(iqs,1,ia) = ysvec(iqs).*phvec(iqs)/vol;
      end
      if (ipref == 1)
         wqmask(iqs,2:4,ia) = 0;
%         wqmask(iqs,2:4,ia) = 0;
%         wqmask(iqs,2:4,ia) = 0;
      end;
      if (idref == 1)
         wqmask(iqs,5:9,ia) = 0;
%         wqmask(iqs,5:9,ia) = 0;
%         wqmask(iqs,5:9,ia) = 0;
%         wqmask(iqs,5:9,ia) = 0;
%         wqmask(iqs,5:9,ia) = 0;
      end;
      %
      iql = find(qvec >=1e-6);
      if (isref == 1)
         wqmask(iql,1,ia) = ysvec(iql).*phvec(iql)/vol;
      end;
      %
      if (ipref == 1)
         wqmask(iql,2,ia) = sqrt(-3)*gkx(inz(iql))./qvec(iql)...
                          .*ypvec(iql).*phvec(iql)/vol;
         wqmask(iql,3,ia) = sqrt(-3)*gky(inz(iql))./qvec(iql)...
                          .*ypvec(iql).*phvec(iql)/vol;
         wqmask(iql,4,ia) = sqrt(-3)*gkz(inz(iql))./qvec(iql)...
                          .*ypvec(iql).*phvec(iql)/vol;
      end;
      %
      if (idref == 1)
         cvec = sqrt(5)*ydvec(iql).*phvec(iql)./qvec(iql).^2/vol;
         wqmask(iql,5,ia) = cvec.*gkx(inz(iql)).*gky(inz(iql))*sqrt(3);
         wqmask(iql,6,ia) = cvec.*gkx(inz(iql)).*gkz(inz(iql))*sqrt(3);
         wqmask(iql,7,ia) = cvec.*gky(inz(iql)).*gkz(inz(iql))*sqrt(3);
         wqmask(iql,8,ia) = cvec.*...
           ((sqrt(3)/2)*(gkx(inz(iql)).^2-gky(inz(iql)).^2));
         wqmask(iql,9,ia) = cvec.*...
           (gkx(inz(iql)).^2+gky(inz(iql)).^2-2*gkz(inz(iql)).^2)/2;
      end;
   end; % for ia = 1:natoms     
   %
   % compress wqmask into a matrix
   %
   wqmat = [];
   for ia = 1:natoms
     for j = 1:9
       if ( sum(abs(wqmask(:,j,ia)))~=0 )
         wqmat = [wqmat wqmask(:,j,ia)]; 
         if (j <= 1)
            % s-state
            wqsign = [wqsign; isNL(1,itype)];
         elseif (j <= 4)
            % p-states
            wqsign = [wqsign; isNL(2,itype)];
         else
            % d-states
            wqsign = [wqsign; isNL(3,itype)];
         end;
       end;
     end;
     %fprintf('ia = %d, size(wqmat,2)=%d\n', ia, size(wqmat,2));
   end;
   %
   % for some reason, I need to conjugate wqmat to get 
   % the same results as what I get from PEtot
   % (mixed up fft convention somewhere?)
   %
   wqmat = conj(wqmat)*sqrt(vol);
end;
fprintf('Time = %11.3e\n', cputime-t0);
