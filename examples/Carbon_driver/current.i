set nrowmax 1;

set cell 12  0  0 0  12   0 0  0 12;

kpoint delete 0 0 0
kpoint add 0.0001 0 0 1.
  
species carbon carbon_pbe.xml

atom C carbon  0.00000000  0.00000000 0.00000000
# atom C carbon 6.00000000 0.00000000 0.00000000

set ecut 35

set xc PBE

set wf_dyn PSD

set ecutprec 20

set nempty 4

set fermi_temp 1000

randomize_wf

load wf.xml

run 0 50 10

save wf.xml

