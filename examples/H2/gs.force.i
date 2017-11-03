set nrowmax 1

#set cell 50 0 0 0 50 0 0 0 50
#set cell 30 0 0 0 30 0 0 0 30
set cell 20 0 0 0 20 0 0 0 20

kpoint delete 0 0 0
kpoint add 0.0001 0.0 0.0 1.
species hydrogen hydrogen_bare_extended_withelectron.xml

atom H1 hydrogen 14.3 15 15
atom H2 hydrogen 15.7 15 15

set atoms_dyn MD

set ecut  30

# should not use ecuts in MPBC
#set ecuts 10

set xc PBE
set wf_dyn PSD

set ecutprec 10

set nempty 2

#randomize_wf

load wf.relaxed.xml

run 0 3 3


