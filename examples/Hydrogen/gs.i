set nrowmax 1

set cell 30 0 0 0 30 0 0 0 30

kpoint delete 0 0 0
kpoint add 0.0001 0.0 0.0 1.
species hydrogen hydrogen_bare_extended.xml

#atom H hydrogen 25 25 25
atom H hydrogen 15 15 15
set ecut  30

# should not use ecuts in MPBC
set ecuts 20

set xc PBE
set wf_dyn PSD

set ecutprec 10

set nempty 2

randomize_wf
save wf_random.xml
load wf.xml

run 0 200 10

save wf.xml

