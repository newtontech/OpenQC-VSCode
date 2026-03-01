%chk=ts_search.chk
%mem=8GB
%nproc=12
# B3LYP/6-31++G(d,p) Opt=(TS,CalcFC) Freq

Transition state search for SN2 reaction

-1 1
C   0.000000  0.000000  0.000000
Cl  1.800000  0.000000  0.000000
H  -0.360000  1.030000  0.000000
H  -0.360000 -0.520000  0.900000
H  -0.360000 -0.520000 -0.900000
F  -2.500000  0.000000  0.000000

