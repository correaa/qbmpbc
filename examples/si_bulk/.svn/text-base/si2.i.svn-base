si2.sys
set ecut 12 
set nempty 8
kp12fcc.i 
# add k-points

# L-Gamma line
kpoint add 0.50  0.50  0.50   0.00
kpoint add 0.45  0.45  0.45   0.00
kpoint add 0.40  0.40  0.40   0.00
kpoint add 0.35  0.35  0.35   0.00
kpoint add 0.30  0.30  0.30   0.00
kpoint add 0.25  0.25  0.25   0.00
kpoint add 0.20  0.20  0.20   0.00
kpoint add 0.15  0.15  0.15   0.00
kpoint add 0.10  0.10  0.10   0.00

# Gamma-X line
kpoint add 0.00  0.00  0.00   0.00
kpoint add 0.05  0.00  0.05   0.00
kpoint add 0.10  0.00  0.10   0.00
kpoint add 0.15  0.00  0.15   0.00
kpoint add 0.20  0.00  0.20   0.00
kpoint add 0.25  0.00  0.25   0.00
kpoint add 0.30  0.00  0.30   0.00
kpoint add 0.35  0.00  0.35   0.00
kpoint add 0.40  0.00  0.40   0.00
kpoint add 0.45  0.00  0.45   0.00
kpoint add 0.50  0.00  0.50   0.00

# X-W line
kpoint add 0.55  0.05  0.50   0.00
kpoint add 0.60  0.10  0.50   0.00
kpoint add 0.65  0.15  0.50   0.00
kpoint add 0.70  0.20  0.50   0.00
kpoint add 0.75  0.25  0.50   0.00

set wf_dyn PSD
set ecutprec 4

randomize_wf
set wf_diag T
set charge_mix_coeff 1.0
run 0 30 5 
save si2.xml 
