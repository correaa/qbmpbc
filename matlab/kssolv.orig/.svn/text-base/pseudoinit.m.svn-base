function pseudovar = pseudoinit(mol)
%
% Usage: pseudovar = pseudoinit(mol);
%
% Purpse:
%    Set up pseudo potential related variables
%    pseudovar contains the following fields:
%
% Input: 
%    mol  --- a Molecule object
% Output:
%    pseudovar --- a structure that contains the following fields: 
%                          qi2
%                          qi3
%                          vq 
%                          rhoq
%                          vqT
%                          rhocq 
%                          Ealphat
%                          wq2 
%

atomlist = get(mol,'atomlist');
pseudovar = [];
if (~isempty(atomlist))
    fprintf('Loading Pseudopotential ...\n');
    t0 = cputime;
    %
    % mnq and mtype are hardcoded in PEToT
    %
    mnq   = 4000;
    mtype = 6;

    pseudovar.qi2     = zeros(mnq,1);
    pseudovar.qi3     = zeros(mnq,1);
    pseudovar.vq      = zeros(mnq,mtype);
    pseudovar.rhoq    = zeros(mnq,mtype);
    pseudovar.vqT     = zeros(mnq,mtype);
    pseudovar.rhocq   = zeros(mnq,mtype);
    pseudovar.Ealphat = 0.0;
    %
    % the number of different atom types
    %
    elements = zeros(110,1);
    ealpha   = zeros(110,1);
    %  
    atomlist = get(mol,'atomlist');
    natoms   = get(mol,'natoms');
    vol      = get(mol,'vol');
    ecut2    = get(mol,'ecut2');
    rcut     = get(mol,'rcut');
    %
    % determine the number of different types of atoms
    %
    for j = 1:natoms
       a = atomlist(j);
       anum = get(a,'anum');
       elements(anum) = 1;
    end;
    ntypes = sum(elements);
    %
    %
    % the dimension of ri and amr, which is 201, are also hardcoded in PEToT  
    %
    ri  = 10*(0:200)/200';
    amr = ones(201,1);
    %
    anums = find(elements==1);
    for j = 1:ntypes
       anum = anums(j);
       atom = Atom(anum);
       fprintf(' atom type = %s\n', get(atom,'symbol'));
       iloc = get(atom,'iloc');
       ic   = get(atom,'ic'); % whether core correction is needed

       if (iloc == 1) 
	 isref = 0;
       elseif (iloc ==2)
	 ipref = 0;
       else
	 idref = 0;
       end;

       icore(j) = ic;
       vwr = loadppt(atom);
       nrr = size(vwr,1);
       rhoc = zeros(nrr,1);

       r    = vwr(:,1);
       vs   = vwr(:,2);
       vp   = vwr(:,3);
       vd   = vwr(:,4);
       ws   = vwr(:,5);
       wp   = vwr(:,6);
       wd   = vwr(:,7);

       if (iloc == 0 & ic == 0)
	  vloc = vwr(:,8);
       elseif (iloc ~= 0 & ic == 1)
	  rhoc = vwr(:,8);
       elseif (iloc == 0 & ic == 1)
	  vloc = vwr(:,8);
	  rhoc = vwr(:,9);
       end;
       %
       if ( iloc == 1)  
	  vloc = vs;
       end;
       if ( iloc == 2) 
	  vloc = vp;
       end;
       if ( iloc == 3)
	  vloc = vd;
       end;
       if ( iloc == 12)
	  vloc = (vs + vp)/2;
       end;
       if ( iloc == 13) 
	  vloc = (vs + vd)/2;
       end;
       if ( iloc == 23) 
	  vloc = (vp + vd)/2;
       end; 

       s = 0.0;
       ch = get(atom,'venum');
       for i = 2:nrr-1
	  if (r(i) < 15.0)
	     s=s+(ch*r(i)+vloc(i)*r(i)^2)*(r(i+1)-r(i-1))/2;
	  end;
       end;
       ealpha(anum) = s*4*pi; 
       %
       % vlocT is for the use of Thomas procedure
       %
       occs = get(atom,'occs');
       occp = get(atom,'occp');
       occd = get(atom,'occd');

       id = find(r<4.0);
       if (~isempty(id)) 
	  a  = vs(id)*occs.*ws(id).^2 + vp(id)*occp.*wp(id).^2 ...
	     + vd(id)*occd.*wd(id).^2;
	  b  = occs*ws(id).^2 + occp*wp(id).^2 + occd*wd(id).^2;
	  vlocT(id) = a./b;
       end;
       id = find(r>=4.0);
       if (~isempty(id))
	  vlocT(id) = vp(id);
       end;
       %
       qmx = sqrt(ecut2)*1.2;
       occt = occs + occp + occd;
       rst = 10.0;   % funny constant;
       %
       %**********************************************************
       %**** for the ionic potential and atomic charge density
       %**** using vs as the local potential for III-V
       %***********************************************************
       %
       iq = 1:mnq;
       g1=(iq-1)*qmx/(mnq-1);
       ismall = find(g1 < 1e-3);
       g1(ismall) = 1e-3; 
       g124pi = 4*pi./g1.^2;

       ifirst = find(r>rst,1);
       idx = 1:ifirst-1;
       ip1 = 2:ifirst;

       r1 = r(idx);
       r2 = r(ip1);
       X1 = r1*g1;
       X2 = r2*g1;
       X21 = X2-X1;
       COSX12 = cos(X1)-cos(X2);
       SINX21 = sin(X2)-sin(X1);

       b1 = vloc(idx).*r1;
       b2 = vloc(ip1).*r2 - b1;

       B1 = repmat(b1,1,mnq);
       B2 = repmat(b2,1,mnq)./X21;
       C1 = COSX12.*B1 - X21.*B2.*cos(X2) + B2.*SINX21;
       s = sum(C1);

       rhoi  = (occs*ws(idx).^2+occp*wp(idx).^2+occd*wd(idx).^2)/(4*pi);
       rhoi1 = (occs*ws(ip1).^2+occp*wp(ip1).^2+occd*wd(ip1).^2)/(4*pi);

       b11 = rhoi.*r1;
       B11 = repmat(b11,1,mnq);
       B22 = repmat(rhoi1.*r2-b11,1,mnq)./X21;
       C11 = COSX12.*B11-X21.*B22.*cos(X2)+B22.*SINX21;
       s1 = sum(C11);

       b1 = vlocT(idx)'.*r1;
       B1 = repmat(b1,1,mnq);
       B2 = repmat(vlocT(ip1)'.*r2-b1,1,mnq)./X21;
       C1 = COSX12.*B1-X21.*B2.*cos(X2)+B2.*SINX21;
       s2 = sum(C1);

       b1 = rhoc(idx).*r1;
       B1 = repmat(b1,1,mnq);
       B2 = repmat(rhoc(ip1).*r2-b1,1,mnq)./X21;
       C1 = COSX12.*B1-X21.*B2.*cos(X2)+B2.*SINX21;
       s3 = sum(C1);

       s2=s2-s;

       s=s.*g124pi;
       s=s-ch*4*pi*cos(g1*r(ifirst))./g1.^2;

       s1=s1.*g124pi;
       s2=s2.*g124pi;
       s3=s3.*g124pi;

       pseudovar.qi2       =g1;
       pseudovar.vq(:,j)   =s;
       pseudovar.rhoq(:,j) =s1;
       pseudovar.vqT(:,j)  =s2;
       pseudovar.rhocq(:,j)=s3;
       clear atom;

       %
       %
       %********************************************************
       %*** for the nonlocal potentail Kleiman-Bylander ref. dv*psi
       %*** using vs as the local potential
       %********************************************************
       %*** isp=1, s state; 2, p state; 3, d state.
       %********************************************************
       rst = rcut/1.03; % funny constant again

       % the following block is taken out of the inner loop
       % begin X setup
       ifirst = find(r>rst,1);
       rv = r(2:ifirst-1);
       rp1 = r(3:ifirst);
       rm1 = r(1:ifirst-2);

       X = rv*g1;
       % end X setup

       for isp = 1:3
	  switch isp
	     case 1
		fprintf('    s-orbit...\n');
	     case 2
		fprintf('    p-orbit...\n');
	     case 3
		fprintf('    d-orbit...\n');
	  end;
	  s = 0.0;
	  s_w = 0.0;
	  for i = 1:nrr
	     if ( r(i) > rst )
		amrI = 0.0;
	     else
		r2 = r(i)/rcut;
		ir = floor(1+r2*200);
		f1 = (ri(ir+1)-r2)/(ri(ir+1)-ri(ir));
		f2 = (r2-ri(ir))/(ri(ir+1)-ri(ir));
		amrI=amr(ir)*f1+amr(ir+1)*f2;
		amrI=1.0/amrI;
	     end;
	     %
	     if (isp == 1)
		vw(i)=(vs(i)-vloc(i))*ws(i);
		vw4(i)=ws(i);
		if (r(i)<rst & i > 1)
		   s=s+ws(i)^2*(vs(i)-vloc(i))*(r(i+1)-r(i-1))/2*r(i)^2;
		   s_w=s_w+ws(i)^2*(r(i+1)-r(i-1))/2*r(i)^2;
		end;
	     end;
	     %
	     if (isp == 2)
		vw(i)=(vp(i)-vloc(i))*wp(i);
		vw4(i)=wp(i);
		if ( (r(i)<rst) & (i>1) ) 
		   s=s+wp(i)^2*(vp(i)-vloc(i))*(r(i+1)-r(i-1))/2*r(i)^2;
		   s_w=s_w+wp(i)^2*(r(i+1)-r(i-1))/2*r(i)^2;
		end;
	     end;
	     %
	     if ( isp == 3)
		vw(i)=(vd(i)-vloc(i))*wd(i);
		vw4(i)=wd(i);
		if ( (r(i)<rst) & (i>1) )
		   s=s+wd(i)^2*(vd(i)-vloc(i))*(r(i+1)-r(i-1))/2*r(i)^2;
		   s_w=s_w+wd(i)^2*(r(i+1)-r(i-1))/2*r(i)^2;
		end;
	     end;
	     vw2(i)=vw(i);
	     vw(i)=vw(i)*amrI;
	  end;
	  %
	  s=4*pi*s;
	  s_w=4*pi*s_w;

	  if( (isp == iloc) | abs(s/s_w) < 1e-10 )
	     scale=0.0;
	  else
	     scale=1/sqrt(abs(s));
	  end;

	  if(s >= 0.0)
	     isNL(isp,j)=1;
	  else
	     isNL(isp,j)=-1;
	  end;
	  %
	  if ( isp == 1 )
              A = sin(X)./X;
%             A = sinc(X/pi);
%             A = besselj(0,X);
	  elseif (isp == 2)
              A = sin(X)./X.^2 - cos(X)./X;
%             A = besselj(1,X);
	  elseif (isp == 3)
              A = (3./X.^3 - 1./X).*sin(X)-3*cos(X)./X.^2;
%             A = besselj(2,X);
	  end;
	  %
	  rv2 = rv.^2.*(rp1-rm1)/2;
	  b = rv2.*vw(2:ifirst-1)';
	  b1 = rv2.*vw2(2:ifirst-1)';
	  b2 = rv2.*vw4(2:ifirst-1)';
	  B = repmat(b, 1, mnq);
	  s = sum(A.*B);
	  B = repmat(b1, 1, mnq);
	  s1 = sum(A.*B);
	  B = repmat(b2, 1, mnq);
	  s2 = sum(A.*B);

	  s = s*4*pi*scale;
	  s1 = s1*4*pi*scale;
	  s2 = s2*4*pi;

	  qi  = g1;
          pseudovar.qi3 = g1;
	  wq(:,isp,j) = s';
	  wq2(:,isp,j) = s1';
	  wq4(:,isp,j) = s2';
       end; % for isp = 1:3
       pseudovar.wq2 = wq2;
       pseudovar.isNL = isNL;
       %
       for i=1:2000
	  rx=(i-1)*10/2000;
	  for i1 = 2:nrr
	     if ( r(i1) > rx ) 
		i1m=i1-1;
		if ( r(i1)-r(i1m) > 1e-10 ) 
		   f1=(r(i1)-rx)/(r(i1)-r(i1m));
		   f2=1.0-f1;
		else
		   f1=0.0;
		   f2=1.0;
		end;
		rho_atomIr(i,j) = ...
		    f1*(occs*ws(i1m)^2+occp*wp(i1m)^2+occd*wd(i1m)^2)/(4*pi) ...
		  + f2*(occs*ws(i1)^2+occp*wp(i1)^2+occd*wd(i1)^2)/(4*pi);
		break;
	     end;
	  end;
       end;
    end; % j = 1:ntypes
    %
    % calculate Ealphat
    %
    ch = 0.0;
    for j = 1:natoms
       atom = atomlist(j);
       anum = get(atom,'anum');
       pseudovar.Ealphat = pseudovar.Ealphat + ealpha(anum);
       ch = ch + get(atom,'venum'); 
    end;
    pseudovar.Ealphat = pseudovar.Ealphat*ch/vol;
    t1 = cputime - t0;
    fprintf('Time = %11.3e\n', t1);
end;
