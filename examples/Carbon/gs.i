# carbon atom

set nrowmax 1

#set cell 30 0 0 0 30 0 0 0 30
#set cell 12 0 0 0 12 0 0 0 12
set cell 16 0 0 0 16 0 0 0 16

kpoint delete 0 0 0
kpoint add 0.0001 0.0 0.0 1.
species carbon carbon_pbe.xml 

#atom C    carbon     0  0  0
#atom C    carbon     6  0  0
#atom C    carbon     0  6  0
#atom C    carbon     0  0  6
atom C    carbon     3  5  7

#set ecut 20
#set ecut 35
set ecut 80

set xc PBE

set wf_dyn PSD

#set ecutprec 20
set ecutprec 70

#set nempty 4
set nempty 2

#set fermi_temp 1000

randomize_wf
save wf_random.xml

#run 0 1 1
#run 0 50 10
run 0 100 10
#run 0 200 10

save wf.xml
