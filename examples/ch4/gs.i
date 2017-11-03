# CH4 Methane Molecule
 
set nrowmax 1

kpoint delete 0 0 0
kpoint add 0.0001 0.0 0.0 1.

ch4.sys
#ch4_rot.sys

set ecut 35

set xc PBE

set wf_dyn PSD

set ecutprec 20

randomize_wf
save wf_random.xml

#run 0 1 1
run 0 50 10

save wf.xml
