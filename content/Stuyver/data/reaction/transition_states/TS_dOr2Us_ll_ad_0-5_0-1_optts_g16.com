%mem=16000MB
%nprocshared=4
# PBE1PBE Def2SVP Freq EmpiricalDispersion=GD3BJ Opt=(TS, CalcFC, NoEigenTest, MaxCycles=100, MaxStep=10, NoTrustUpdate, RecalcFC=30) Geom=ModRedun NoSymm scrf=(smd,solvent=Water)

 TS_dOr2Us_ll_ad_0-5_0-1_optts_g16

-1 1
C    4.60610140  -0.02730503   0.00493278 
Cl   6.79848937   0.07640948  -0.00821413 
H    4.54314962  -0.57788310  -0.91000959 
H    4.55684153  -0.54992160   0.93698707 
H    4.47363589   1.03400686  -0.01029273 
F    2.57786220  -0.12371660   0.02996660 

B 1 2
B 1 6



