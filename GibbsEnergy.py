T = 100;

TKelvin = T + 273.15;

# Fe
# CP (T) = 28.18 - 7.32e-3 *T -0.29e-6/T^2 + 25.04e6 *T^2
CP0298  = 1/(TKelvin - 298) * (28.18 * (TKelvin- 298) + (-7.32e-3) / 2 * (TKelvin - 298)* (TKelvin - 298) + (-0.29e-6) * (1/298 - 1/TKelvin) + 25.04e6/3 * ( 1/TKelvin^3 -  1/(298^3))) # J/mol.K;

# Fe3O4
# CP (T) = 2659.108-2521.53e-7*T + 20.7347e-6*T^2 + 1.36769e-3/T^2 -3.645541e4* T^(-.5) J/mol.K 
CP0Fe3O4 =  1/(TKelvin- 298)*(2659.108 * (TKelvin - 298) - 2521.53e-7/2 *(TKelvin^2 - 298^2) + 20.7347e-6/3 *(TKelvin^3 - 298^3) +  (-3.645541e4) * 2  (math.sqr(TKelvin) - math.sqr(298)))

# H2O
# J/mol.K; Used the value of Cp at 25 C, small change over the temperature of 25-300 C
CP0H2O = 75.2 

# H2 (g)
# (6.52+.78e-3*TKelvin+.12e5/TKelvin^2)* 4.18  can be used for Cp but variation is very small
# J/mol.K; Used the value of Cp at 25 C, small change over the temperature of 25-300 C
CP0H2 = 29.0 


# H+
# J/mol.K; Used the value of Cp at 25 C
CP0Hp = -71.0


#Fe(OH)2
# J/mol.K; Used the averaged value of Cp from Tremaine 
CP0FeOH2 = 133.0 


