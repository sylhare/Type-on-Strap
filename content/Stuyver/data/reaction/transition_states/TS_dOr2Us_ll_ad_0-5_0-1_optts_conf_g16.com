%mem=16000MB
%nprocshared=4
# PBE1PBE Def2SVP Freq EmpiricalDispersion=GD3BJ Opt=(TS, CalcFC, NoEigenTest, MaxCycles=100, MaxStep=10, NoTrustUpdate, RecalcFC=30) Geom=ModRedun NoSymm scrf=(smd,solvent=Water)

 TS_dOr2Us_ll_ad_0-5_0-1_optts_conf_g16

-1 1
C    5.39022900   0.66593400  -1.05053800 
Cl   6.99843600  -0.38338800  -1.97717900 
H    5.05843700   1.16140900  -1.95593900 
H    4.86964700  -0.22373600  -0.71422900 
H    6.01471100   1.20890000  -0.35017000 
F    3.88514000   1.64894100  -0.18283500 

B 1 2
B 1 6



