set nrowmax 1

#set cell 50 0 0 0 50 0 0 0 50
set cell 30 0 0 0 30 0 0 0 30

kpoint delete 0 0 0
kpoint add 0.0001 0.0 0.0 1.
species hydrogen hydrogen_bare_extended.xml

atom H1 hydrogen 14.3 15 15
atom H2 hydrogen 15.7 15 15
set ecut  20

# should not use ecuts in MPBC
#set ecuts 10

set xc PBE
set wf_dyn PSD

set ecutprec 10

set nempty 2

randomize_wf
#load wf.xml #if available

run 0 100 10
save wf.xml

#run 0 100 10
#save wf.xml

