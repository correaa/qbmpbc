set cell 3 0 0 0 50 0 0 0 50
kpoint delete 0 0 0
kpoint add 0.0001 0.0 0.0 1.
set ecut 2

# should not set ecuts in MPBC
#set ecuts 1

set xc PBE
set wf_dyn PSD
#set ecutprec 0.4
set ecutprec 1.0

set nempty 6
randomize_wf
#load wf.xml

run 0 50 10

save wf.xml


