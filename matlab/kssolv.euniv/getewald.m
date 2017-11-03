function Ewald = getewald(mol)
%
% Usage: Ewald = getewald(mol);
%
% Purpose: 
%    The function returns the Ewald energy.
%
% Input:
%    mol   --- a Molecule object.
%
% Output:
%    Ewald --- Ewald energy
%

% check the validity of the input argument.
%
fprintf('Calculating Ewald energy...\n');
atomlist = get(mol,'atomlist');
t0       = cputime;

if (~isempty(atomlist))
    natoms   = get(mol,'natoms');
    xyzlist  = get(mol,'xyzlist');
    %

    C  = get(mol,'supercell');
    vol = get(mol,'vol');
    %
    d1=norm(C(:,1));
    d2=norm(C(:,2));
    d3=norm(C(:,3));
    %
    dave=(d1*d2*d3)^(1/3);
    if (dave <= 16)
       beta = 3.0/dave;
    else
       beta = 5.0/dave;
    end;
    %
    nc1=floor(16.0/beta/d1+1);
    nc2=floor(16.0/beta/d2+1);
    nc3=floor(16.0/beta/d3+1);
    %
    % retrieve reciprocal space info, use ecut2
    %
    ecut2 = get(mol,'ecut2');
    gmask = FreqMask(mol,ecut2);
    gkk2 = get(gmask,'gkk');
    gm2  = get(gmask,'gm');
    gkx2 = get(gmask,'gkx');
    gky2 = get(gmask,'gky');
    gkz2 = get(gmask,'gkz');
    %
    [X,Y,Z]=meshgrid(-nc1:nc1,-nc2:nc2,-nc3:nc3);
    [n1,n2,n3]=size(X);
    ijkmat = [reshape(X,1,n1*n2*n3); ...
              reshape(Y,1,n1*n2*n3); ...
              reshape(Z,1,n1*n2*n3)];
    %
    for ia = 1:natoms
       atomchg(ia) = get(atomlist(ia),'venum');
    end;
    %
    Ewald = 0.0;
    for ia1 = 1:natoms
       fatom(1,ia1) = 0.0;
       fatom(2,ia1) = 0.0;
       fatom(3,ia1) = 0.0;

       xyz1 = xyzlist(ia1,:)';

       for ia2 = 1:natoms

	  xyz2 = xyzlist(ia2,:)';

	  dxyz0 = xyz2 - xyz1;

	  ch1=atomchg(ia1);
	  ch2=atomchg(ia2);

	  sr =0.0;
	  srxyz = zeros(3,1);

%         for i = -nc1:nc1
%            for j = -nc2:nc2
%               for k = -nc3:nc3
%                  xyz = dxyz0 + C*[i j k]';
%                  r=norm(xyz);
%
%		   if (beta*r > 8.0) 
%                     continue;
%                  end;
%                  if (ia2 ~= ia1) 
%                     derfrdr=(erfc(beta*(r+2.0e-5))/(r+2.0e-5)...
%                            -erfc(beta*r/r))/2.0e-5;
%
%                     srxyz = srxyz - derfrdr*xyz/r;
%                  end; 
%                  if (ia1 == ia2 & i == 0 & j == 0 & k == 0) 
%                     continue;
%                  end; 
%                  sr=sr+erfc(beta*r)/r;
%               end; % k  
%            end; % j 
%         end; % i

          xyzmat = repmat(dxyz0,1,n1*n2*n3) + C*ijkmat;
          rvec = sqrt(sum(xyzmat.*xyzmat,1));
          iz = find(rvec==0);
          if (~isempty(iz))
             rvec(iz) = [];
          end;
          i8 = find(rvec<=8.0/beta);
          if ( ia1 ~= ia2)
             vderfrdr=(erfc(beta*(rvec(i8)+2.0e-5))./(rvec(i8)+2.0e-5)...
                           -erfc(beta*rvec(i8)./rvec(i8)))/2.0e-5;
             srxyzmat = - (ones(3,1)*vderfrdr).*xyzmat(:,i8)...
                          ./(ones(3,1)*rvec(i8));
             srxyz = sum(srxyzmat,2);
          end
          sr = sum(erfc(beta*rvec(i8))./rvec(i8));
	  if (ia1 == ia2)
	     sr=sr-2*beta/sqrt(pi);
	  end;
%
%         factor 2 for the sum 1,ng2, is only half the sphere
%
	  Ewald=Ewald+sr*ch1*ch2*0.5;

	  fatom(1,ia1)=fatom(1,ia1)+srxyz(1)*ch1*ch2;
	  fatom(2,ia1)=fatom(2,ia1)+srxyz(2)*ch1*ch2;
	  fatom(3,ia1)=fatom(3,ia1)+srxyz(3)*ch1*ch2;
       end; % ia2 
    end; % ia1
    %
    sq  = 0.0;
    inz = find( (abs(gkk2) > 1e-10) & (abs(gkk2) <= 50) );
    nnz = length(inz);
    gkk2 = gkk2*2;  % scale by a factor of 2 to be consistent with PEtot
    yy = gkk2(inz)/(4*beta^2);
    ff1 = exp(-yy);
    ff2 = (ff1*4*pi/vol*2)./gkk2(inz);
%    for i = 1:nnz
%       scos = 0.0;
%       for ia1 = 1:natoms
%	  ssin = 0.0; 
%	  xyz1 = xyzlist(ia1,:)';
%
%	  for ia2 = 1:natoms
%	     xyz2 = xyzlist(ia2,:)';
%	     dxyz0 = xyz2-xyz1; 
%
%	     ph=[gkx2(inz(i)) gky2(inz(i)) gkz2(inz(i))]*dxyz0;
%
%	     ch=atomchg(ia1)*atomchg(ia2);
%
%	     scos=scos+ch*cos(ph);
%	     ssin=ssin+ch*sin(ph);
%	  end; % ia2
%
%	  ssin=ssin*ff2(i);
%	  fatom(1,ia1)=fatom(1,ia1)+ssin*gkx2(inz(i));
%	  fatom(2,ia1)=fatom(2,ia1)+ssin*gky2(inz(i));
%	  fatom(3,ia1)=fatom(3,ia1)+ssin*gkz2(inz(i)); 
%       end; % ia1
%       sq = sq + ff1(i)*scos/gkk2(inz(i));
%    end;
%
%   ----- The OUTER LOOP (for i=1:nnz) above has been vectorized ---
    chcosvec = 0.0;
    
    gkxyz2 = zeros(length(inz),3);
    gkxyz2(:,1) = gkx2(inz);
    gkxyz2(:,2) = gky2(inz);
    gkxyz2(:,3) = gkz2(inz);

    for ia1 = 1:natoms
       chsinvec = 0.0;
       xyz1  = xyzlist(ia1,:)';
       for ia2 = 1:natoms
          xyz2  = xyzlist(ia2,:)';
          dxyz0 = xyz2 - xyz1;
          phvec = gkxyz2*dxyz0;
          ch=atomchg(ia1)*atomchg(ia2);
          chcosvec = chcosvec + ch*cos(phvec);
          chsinvec = chsinvec + ch*sin(phvec);
       end;
       fssin = chsinvec.*ff2;
       fatom(1,ia1) = fatom(1,ia1)+fssin'*gkx2(inz);
       fatom(2,ia1) = fatom(2,ia1)+fssin'*gky2(inz);
       fatom(3,ia1) = fatom(3,ia1)+fssin'*gkz2(inz);
    end 
    sq = sum(ff1.*chcosvec./gkk2(inz));

    sq = sq*4*pi/vol; % don't need to divide 2 because we've used the whole sphere
    ch = sum(sum(atomchg'*atomchg));
    sq = sq - ch*pi/beta^2/vol;
    Ewald = Ewald + sq/2;
else
    Ewald = 0.0; 
end;
fprintf('Time = %11.3e\n', cputime - t0);
