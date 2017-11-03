% plot NMR shift convergence with respect to ecut
%  sigma (ppm) of H2 molecule (oriented parallel/perpendicular to B)

data = [ % from sigma.dat
%    X_old    X        X_fast     Y        Y_fast
 30  25.973   25.7785  25.7789    23.1863  23.1929
 40  26.251   25.959   25.9591    23.6386  23.6444 
%50  27.266   26.5525  26.8063    30.9408  31.3387 
 60  26.684   26.6712  26.6713    24.8338  24.8361 
 70  26.922   26.8968  26.8967    25.2066  25.2087 
 80  27.051   0        27.0276    25.3845  25.3854 
 90  27.117   0        27.0922    25.3761  25.3748  
%100 27.703   27.3774  27.4784    0        27.6208
%110 0        0        27.4781    0        27.6459
 140 0        27.4328  27.4333    0        25.5472
 200 0        0        27.5119    0        25.7483
 260 0        0        27.5674    0        25.8296
%320 0        0        0          0        0
];

Emax = 350;
fs = 14;
        
eng_data = [
%ecut  E_X_fast           E_Y_fast
 30    -2.59160130314592 -2.46677337472247
 40    -2.66246390043303 -2.53269101526979
 60    -2.76892338215787 -2.63529423397661
 70    -2.79044065136564 -2.65614866017335
 80    -2.80985967503119 -2.67488438821231
 90    -2.84211240025006 -2.70566371532681
 140   -2.89179680657282 -2.75136000629495
 200   -2.90442692719796 -2.76274414759549
 260   -2.91191140716424 -2.77009342777496
];

sig_para_anl = 27.17; % ppm
sig_perp_anl = 25.77; % ppm
sig_avg_anl  = (sig_para_anl + sig_perp_anl*2)/3;

figure(3); 
plot(data(:,1), data(:,4)-data(:,6), 'b.-', ...
     [0 Emax],(sig_para_anl-sig_perp_anl)*[1 1], 'b--'); 
set(gca,'FontName','Times','FontSize',fs);
xlabel('E_{cut}   (Ry)');
ylabel('\sigma_{||} - \sigma_{\perp}   (ppm)');
title('NMR shift of H_2 molecule');
xlim([0 Emax]); 
ylim([1  3])

figure(4); 
plot(data(:,1), data(:,4), 'b.-', ...
     data(:,1), data(:,6), 'm+-', ...
     data(:,1), (data(:,4)+data(:,6)*2)/3, 'ks-', ...
     [0 Emax],sig_para_anl*[1 1], 'b--', ...
     [0 Emax],sig_perp_anl*[1 1], 'm--', ...
     [0 Emax],sig_avg_anl*[1 1], 'k--' ); 
set(gca,'FontName','Times','FontSize',fs);
%text(105,26.8,'27.17','FontSize',fs-1);
xlabel('E_{cut}   (Ry)');
ylabel('\sigma   (ppm)');
legend('\sigma_{||}','\sigma_{\perp}','\sigma_{avg}',4);
title('NMR shift of H_2 molecule');
xlim([0 Emax]); 
ylim([23 28])
%print -depsc H2_NMR_converge

figure(5);
subplot(2,1,1);
plot(eng_data(:,1), eng_data(:,2), 'b.-', ...
     eng_data(:,1), eng_data(:,3), 'm+-' );
set(gca,'FontName','Times','FontSize',fs);
legend('E_{||}', 'E_\perp');
title('Energy of H_2 molecule in magnetic field');
ylabel('E (Hartree)');
xlim([0 Emax]);
subplot(2,1,2);
plot(eng_data(:,1), eng_data(:,3)-eng_data(:,2), 'k.-');
set(gca,'FontName','Times','FontSize',fs);
%legend('\Delta E',4);
xlabel('E_{cut}   (Ry)');
ylabel('\Delta E (Hartree)');
xlim([0 Emax]); 
%print -depsc H2_eng_converge

sig_avg_num = (data(:,4)+data(:,6)*2)/3

fprintf('sigma (anl) = %.2f  sigma (num) = %.2f  error = %e\n', ...
    sig_avg_anl, sig_avg_num(end), (sig_avg_num(end)-sig_avg_anl)/sig_avg_anl);

