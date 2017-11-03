# compute the stress tensor for FCC Au
set nrowmax 16
# FCC Au exp lattice constant a=4.08 Ang = 7.71 Bohr
set cell  3.9 3.9 0.000  0.000 3.9 3.9  3.9 0.000 3.9
set ref_cell 4 4 0  0 4 4  4 0 4
species gold au-pbe.TM.xml
atom Au gold 0 0 0

set nempty 8
set ecut 90

# define kpoints using file kp8881fcc.i
kp8881fcc.i

set stress ON
set ecuts 85
set xc PBE
set wf_dyn PSD
set ecutprec 10
randomize_wf
set charge_mix_coeff 0.2
set fermi_temp 300
run 10 5 5
# save sample for later use
save gs.xml
