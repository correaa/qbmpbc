function [vion,vionT,rho] = getvion(mol,pseudovar)
%
% Usage: [vion,vionT,rho] = getvion(mol,pseudovar);
%
% Purpose:
%    Computes the local ionic potential from ionic pseudopotential
%
% Input:
%    mol       --- a Molecule object
%    pseudovar --- pseudopotential information
% Output:
%    vion      --- 3D local ionic potential
%    vionT     --- currently not used
%    rho       --- 3D charge density associated with vion
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% have to check the input argument

fprintf('Calculating local ionic potential...\n');
t0       = cputime;
n1       = get(mol,'n1');
n2       = get(mol,'n2');
n3       = get(mol,'n3');
vol      = get(mol,'vol');
atomlist = get(mol,'atomlist');

vion  = zeros(n1,n2,n3);
vionT = zeros(n1,n2,n3);
rho   = zeros(n1,n2,n3);

if (~isempty(atomlist))
    xyzlist  = get(mol,'xyzlist');
    natoms   = get(mol,'natoms');
    anumlist = zeros(natoms,1);


    vionf  = zeros(n1,n2,n3);
    vionTf = zeros(n1,n2,n3);
    rhof   = zeros(n1,n2,n3);
    rhocr  = zeros(n1,n2,n3);
    n123   = n1*n2*n3;

    qi2   = pseudovar.qi2';
    vq    = pseudovar.vq;
    vqT   = pseudovar.vqT;
    rhoq  = pseudovar.rhoq;
    rhocq = pseudovar.rhocq;

    totnel = 0;
    elements = zeros(110,1);

    for j = 1:natoms
       a = atomlist(j);
       anum = get(a,'anum');
       elements(anum) = 1;
       %
       % count the total number of valence electrons
       %
       totnel = totnel + get(a,'venum');
       %
       % precompute atomic number list 
       %
       atom = atomlist(j);
       anumlist(j) = get(atom,'anum');
    end;
    ntypes = sum(elements);
    ic     = zeros(ntypes,1);
    %
    icc = 0;
    %
    % list of atomic numbers (no repeat)
    %
    anums = find(elements==1);
    for itype = 1:ntypes
       atom = Atom(anums(itype));
       icc = icc + get(atom,'ic');
       ic(itype) = get(Atom(anums(itype)),'ic');
    end;
    %
    ecut2  = get(mol,'ecut2');
    gmask2 = FreqMask(mol,ecut2);
    ng2    = get(gmask2,'ng');
    gm2    = get(gmask2,'gm');
    gkk2   = get(gmask2,'gkk');
    gkx2   = get(gmask2,'gkx');
    gky2   = get(gmask2,'gky');
    gkz2   = get(gmask2,'gkz');
    inz    = find(gm2~=0);

    mnq = 4000;  % a constant set in param.escan

    % vectorize the original for i = 1:nnz loop
    for itype = 1:ntypes
       index  = find(anumlist == anums(itype));
       xyzmat = xyzlist(index,:)';
       phmat  = [gkx2(inz) gky2(inz) gkz2(inz)]*xyzmat;
       ccvec  = sum(exp(sqrt(-1).*phmat),2);
       %
       qvec  = sqrt(gkk2(inz)*2);
       iqvec = floor(1+qvec*(mnq-1)/qi2(mnq));
       %
       xvec  = (qvec-qi2(iqvec))./(qi2(iqvec+1)-qi2(iqvec));
       %
       % some kind of interpolation) ....
       %
       f1vec = 1-xvec-xvec.*(1-xvec)/2;
       f2vec = xvec + xvec.*(1-xvec);
       f3vec = -xvec.*(1-xvec)/2;
       %
       yvec = vq(iqvec,itype).*qi2(iqvec).^2.*f1vec   ...
	   + vq(iqvec+1,itype).*qi2(iqvec+1).^2.*f2vec ...
	   + vq(iqvec+2,itype).*qi2(iqvec+2).^2.*f3vec;

       y2vec = vqT(iqvec,itype).*qi2(iqvec).^2.*f1vec ...
	     + vqT(iqvec+1,itype).*qi2(iqvec+1).^2.*f2vec ...
 	     + vqT(iqvec+2,itype).*qi2(iqvec+2).^2.*f3vec;
       %
       iqs = find(qvec<1e-6);
       yvec(iqs) = 0;
       y2vec(iqs) = vqT(4,itype);  
       %
       iql = find(qvec>=1e-6);
       yvec(iql) = yvec(iql)./qvec(iql).^2;
       y2vec(iql) = y2vec(iql)./qvec(iql).^2;
       %
       y1vec = rhoq(iqvec,itype).*f1vec+rhoq(iqvec+1,itype).*f2vec+rhoq(iqvec+2,itype).*f3vec;
       %
       vionf(inz)  = vionf(inz) + yvec.*ccvec/vol;
       vionTf(inz) = vionTf(inz)+ y2vec.*ccvec/vol;
       rhof(inz)   = rhof(inz)  + y1vec.*ccvec/vol;
       %
       if ( ic(itype) ~= 0 )
	  y1vec = rhocq(iqvec,itype).*f1vec+rhocq(iqvec+1,itype).*f2vec+rhocq(iqvec+2,itype).*f3vec;
	  rhocrf= rhocrf+y1vec.*ccvec/vol;
       end
    end;

%    nnz = length(inz);
%    for i = 1:nnz
%      for itype = 1:ntypes
%	 %cc = 0;
%	 %for ia = 1:natoms
%	 %  atom = atomlist(ia);
%	 %  anum = get(atom,'anum');
%	 %   if ( anum == anums(itype) )
%	 %      xyzr = xyzlist(ia,:)';
%	 %      ph   = [gkx2(inz(i)) gky2(inz(i)) gkz2(inz(i))]*xyzr;
%	 %      cc   = cc+exp(sqrt(-1)*ph);
%	 %   end;
%	 %end;
%	 % vector version
%	 index = find(anumlist == anums(itype));
%	 xyzr = xyzlist(index,:)';
%	 ph   = [gkx2(inz(i)) gky2(inz(i)) gkz2(inz(i))]*xyzr;
%	 cc = sum(exp(sqrt(-1).*ph));     
%	 %
%	 q  = sqrt(gkk2(inz(i))*2);
%	 iq = floor(1+q*(mnq-1)/qi2(mnq));
%	 %
%	 x  = (q-qi2(iq))/(qi2(iq+1)-qi2(iq));
%	 %
%	 % some kind of interpolation) ....
%	 %
%	 f1 = 1-x-x*(1-x)/2;
%	 f2 = x + x*(1-x);
%	 f3 = -x*(1-x)/2;
%	 %
%	 y = vq(iq,itype)*qi2(iq)^2*f1   ...
%	   + vq(iq+1,itype)*qi2(iq+1)^2*f2 ...
%	   + vq(iq+2,itype)*qi2(iq+2)^2*f3;
%
%	 y2 = vqT(iq,itype)*qi2(iq)^2*f1 ...
%	    + vqT(iq+1,itype)*qi2(iq+1)^2*f2 ...
%	    + vqT(iq+2,itype)*qi2(iq+2)^2*f3;
%
%	 if (q < 1e-6)
%	    y  = 0;
%	    y2 = vqT(4,itype);
%	 else
%	    y2 = y2/q^2;
%	    y  = y/q^2; 
%	 end;
%	 %
%	 y1 = rhoq(iq,itype)*f1+rhoq(iq+1,itype)*f2+rhoq(iq+2,itype)*f3;
%
%	 vionf(inz(i))  = vionf(inz(i)) + y*cc/vol;
%	 vionTf(inz(i)) = vionTf(inz(i))+ y2*cc/vol;
%	 rhof(inz(i))   = rhof(inz(i))  + y1*cc/vol;
%    %
%    %    JCM: Not even sure that these values are used, but I'll keep them in here for now 
%    %
%	 if ( ic(itype) ~= 0 )
%	    y1 = rhocq(iq,itype)*f1+rhocq(iq+1,itype)*f2+rhocq(iq+2,itype)*f3;
%	    rhocrf(i)=rhocrf(i)+y1*cc/vol;
%	 end
%       end % for itype = 1:ntypes
%    end; % for i = 1:nnz
%    toc
    %
    % get rho in real space (use fftn instead of ifftn?)
    %
    vion  = real(fftn(vionf));
    rho   = real(fftn(rhof));
    vionT = real(fftn(vionTf));
    %
    % vionT is still not correct!!!!!
    %
    %reshape(vionT,n123,1)
    %pause

    if (icc > 0)
	%  work on rhocr
    end;

    s = sum3d(rho);
    s = s*vol/n123;

    %
    % total number of electron is hard coded now.
    %
    s = totnel/s;
    rho = rho*s;
end;
fprintf('Time = %11.3e\n', cputime-t0);
