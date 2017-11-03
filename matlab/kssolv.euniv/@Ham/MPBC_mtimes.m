function [HX] = MPBC_mtimes(H,X)

if (nargin == 2)

    n1 = get(X,'n1');
    n2 = get(X,'n2');
    n3 = get(X,'n3');
    psi = get(X,'psi');
    idxnz = get(X,'idxnz');
    ncols = get(X,'ncols');
    C = get(H,'supercell');
    CI = inv(C);
    vtot = get(H,'vtot');
    [m1,m2,m3] = size(vtot);
    %vtot = zeros([m1,m2,m3]);  % for Landay levels
    wqmat = [];
    wqsign = [];
    if (m1==n1 & m2==n2 & m3==n3)
        % define the sapce variables
        gXhat = 2*pi*CI(1,1)*[0:n1*n2/2-1, -n1*n2/2:-1]/n2;
        gY = 2*pi*CI(2,2)*[0:n2/2-1,-n2/2:-1];
        if n3 == 1
            gZ = 0;
        else
            gZ = 2*pi*CI(3,3)*[0:n3/2-1,-n3/2:-1];
        end
        gx = C(1,1)*[0:n1-1]/n1;
        gxhat = zeros(1,n1*n2);
        l=1;
        for j=1:n2,
            for i=1:n1,
                gxhat(l) = 2*pi/C(1,1)/C(2,2)*gx(i) + gY(j);
                l = l + 1;
            end
        end
        
        % compute energy
        for m=1:ncols
            psi_XhatZ = zeros(n1*n2,n3);
            psi_XhatZ(idxnz) = psi{m};
            
            % transform from Reciprocal space to Intermediate space
            psi_xhatZ = ifft(psi_XhatZ,[],1);
            
            % unfold in Intermediate space
            psi_xYZ = reshape(psi_xhatZ,n1,n2,n3);
            
            % transform from Intermediate space to Real space
            psi_xyz = ifft(ifft(psi_xYZ,[],2),[],3);

            % apply potential energy in real space
            Vpsi_xyz = vtot.*psi_xyz;
            
            % transform from Real space to Intermediate space
            Vpsi_xYZ = fft(fft(Vpsi_xyz,[],2),[],3);
            
            % fold in Intermediate space
            Vpsi_xhatZ = reshape(Vpsi_xYZ,n1*n2,n3);
            
            % transform from Intermediate space to Reciprocal space
            Vpsi_XhatZ = fft(Vpsi_xhatZ,[],1);
            
            % kinetic energy terms
            T1    = (gXhat'*ones(1,n3)).^2.*psi_XhatZ/2; 
            T2    = fft( (gxhat'*ones(1,n3)).^2.*psi_xhatZ/2,[],1 );
            T3    = (ones(n1*n2,1)*gZ).^2.*psi_XhatZ/2;
            Kpsi_XhatZ = T1 + T2 + T3;

            % Kinetic energy + Potential energy
            HXpsi{m} =Kpsi_XhatZ(idxnz) +  Vpsi_XhatZ(idxnz);
            
            if 0
                if (m==1)
                    figure(1)
                    plot(abs(psi_XhatZ))
                    figure(2)
                    mesh(gXhat'*ones(1,n3))
                    figure(3)
                    plot(abs(psi_xhatZ))
                    figure(4)
                    mesh(gxhat'*ones(1,n3))
                    drawnow
                end
            end
            %
            % apply nonlocal pseudopotential
            wqmat  = get(H,'wqmat');
            wqsign = get(H,'wqsign');
            if (~isempty(wqmat))
                HXpsi{m} = HXpsi{m} + wqmat*(wqsign.*(wqmat'*psi{m}));
            end;
        end
        HX = Wavefun(HXpsi,n1,n2,n3,idxnz);
        %disp('HX information');
        %HX
    else
        error('The dimension of the Hamiltonian does not match that of the wave function');
    end
else
    error(sprintf('Ham: multiplication syntax: Y=H*X'))
end;
    
    









