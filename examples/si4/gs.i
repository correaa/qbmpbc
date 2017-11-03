# ground state of Si4
set cell 20 0 0 0 20 0 0 0 20
species silicon silicon.xml
atom Si1 silicon 3.500 0.000 0.000
atom Si2 silicon 0.000 2.000 0.000
atom Si3 silicon -3.500 0.000 0.000
atom Si4 silicon 0.000 -2.000 0.000
set ecut 12
set wf_dyn PSDA
set ecutprec 3
randomize_wf
run 0 200  
save si4gs.xml
