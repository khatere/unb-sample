import numpy as np
import matplotlib.pyplot as plt
import KeqatT
import math
from scipy.optimize import fsolve
deg = 5                      # Order of the fitted polynomial
T= 25                      # The temperature at which the K is needed to be calculated
Kmorpholine = KeqatT.logKFitT(deg,T)
Kw = KeqatT.KwatT(T)
KFeOH3 = KeqatT.KFeOH3atT(T)
KFeOH2 = KeqatT.KFeOH2atT(T)
KFeOH = KeqatT.KFeOHatT(T)
DebyeHuckConst = KeqatT.DebyeHuckelConst(T)
print  Kmorpholine, Kw, DebyeHuckConst, KFeOH3, KFeOH2, KFeOH

#30-60 mg/kg of added Morpholine (morpholine molar mass: 87.1 g/mol) equals to  0.000344-0.000688 mol/L take an average value of 0.000516 mol/L
ConcC4H9ONTotal = 0.000516 #mol/L

#solve the following equations for ConcH, ConcOH, ConcC4H9ONH, IStrength, gamma1
# Kw = ConcH*gamma1*ConcOH*gamma1;   
# ConH + ConcC4H9ONH = ConOH;  charge neutrality  
# Kmorpholine = ConcC4H9ONH*gamma1*ConOH*gamma1/ConcC4H9ON;
# ConcC4H9ONTotal = ConcC4H9ONH + ConcC4H9ON = ConcC4H9ONH + gamma1^2*ConcC4H9ONH*ConOH/Kmorpholine
# IStrength = ((z^2)*ConcH + (z^2)*ConcOH +(z^2)*ConcC4H9ONH)/2; z=1;
# log(gamma1) =  -DebyeHuckconst*z^2*[(sqrt(I)/(1+sqrt(I)))-Beta*I];   (Davis equation)         
#Since the concentrations are very small, I multiply them by 1e5 and introduced new concentrations and IStrength.
def equations(p):
    out = [(p[0]/100000) - (p[1]/100000) + (p[2]/100000)]
    out.append(((p[0]/100000) + (p[1]/100000) + (p[2]/100000))/2 - (p[3]/100000))
    out.append((p[0]/100000)*(p[1]/100000)*p[4]*p[4] - Kw)
    out.append((p[2]/100000) + (p[2]/100000)*(p[1]/100000)*p[4]*p[4]/Kmorpholine - ConcC4H9ONTotal)
    out.append(math.pow(10,(-1*DebyeHuckConst*(math.sqrt(p[3]/100000)/(1 + math.sqrt(p[3]/100000)) - 0.2*p[3]/100000))) - p[4])
    return out
P1 =  fsolve(equations, [10, 10, 10, 10, 10]) # fsolve is sensitive to initial guess
[ConcH, ConcOH, ConcC4H9ONH, IStrength, gamma1] = [ P1[0]*1e-5, P1[1]*1e-5, P1[2]*1e-5, P1[3], P1[4]] 
print 'ConcH=', ConcH, 'ConcOH=', ConcOH, 'ConcC4H9ONH=', ConcC4H9ONH, 'IStrength=', IStrength, 'gamma1=', gamma1
pH = -1*math.log(ConcH,10)
print 'pH=', pH

# At the oxide-Bulk interface O_B

n = 2;
F = 96485;                   # C/mol
R = 8.314;                   # (j/mol.K)
TKelvin = T + 273.15; 
h = 6.62606957e-34           # Plank's constant J.s
kboltzman = 1.3806488e-23;   # Boltzman constant J/K
Beta =0.5;

Ce_OB # iron species concentration at the O_B interface
CSat_OB


# Corrosion reaction at the M_P interface; All the reaction are considered in the following direction: R -> O + ne-
# Fe + 2H2O = 2H+ +2e- + Fe(OH)2

EeFe_Fe2p_MP = -1*DeltaGFe_Fe2_MP/(n*F) - math.log(10)*2*R*TKelvin/(n*F)*pH + R*TKelvin/(n*F)*math.log(Ce_MP); # Equilibrium potential for the oxidation reaction at M_P interface
i0Fe_Fe2p_MP = F*kboltzman *TKelvin/h*(math.exp(-1*DeltaGFe_Fe2_MP/(R*TKelvin))*Ce_MP*math.exp(Beta*n*F*EeFe_Fe2p_MP/(R*TKelvin)));
iFe_Fe2p_MP = i0Fe_Fe2p_MP *(math.exp(Beta*n*F/(R*TKelvin)*(Ecorr-EeFe_Fe2p_MP)) - math.exp(-1*(1-Beta)*n*F/(R*TKelvin)*(Ecorr-EeFe_Fe2p_MP)));

#Precipitation of magnetite at M_P interface  3Fe(OH)2 = Fe3O4 + 2H2O + 2H+ + 2e-

EeFe2p_Fe3O4_MP = -1*DeltaGFe2p_Fe3O4_MP/(n*F) - math.log(10)*2*R*TKelvin/(n*F)*pH - 3*R*TKelvin/(n*F)*math.log(Ce_MP);

#Hydrogen production 2H+ + 2e- = H2 at M_P; The potential and current equations are written for the R => O + ne- ( H2 = 2H+ + 2e-)
#Standard DeltaGH2_Hp is zero at any temperature

EeH2_Hp_MP =  - math.log(10)*2*R*TKelvin/(n*F)*pH - R*TKelvin/(n*F)*math.log(CeH2_MP); # At M_P interface
i0H2_Hp_MP = F*kboltzman *TKelvin/h*CeH2_MP*math.exp(Beta*n*F*EeH2_Hp_MP/(R*TKelvin));
iH2_Hp_MP = i0H2_Hp_MP *(math.exp(Beta*n*F/(R*TKelvin)*(Ecorr-EeH2_Hp_MP)) - math.exp(-1*(1-Beta)*n*F/(R*TKelvin)*(Ecorr-EeH2_Hp_MP)));

# Mixed Potential at the M_P interface
Ecorr = R*Tkelvin/(Beta-(1-Beta)*n*F)*math.log((i0H2_Hp_MP*math.exp(Beta*n*F*EeH2_Hp_MP/(R*TKelvin))+i0Fe_Fe2p_MP*math.exp(Beta*n*F*EeFe_Fe2p_MP/(R*TKelvin)))/(i0H2_Hp_MP*math.exp(-1*(1-Beta)*n*F*EeH2_Hp_MP/(R*TKelvin))+i0Fe_Fe2p_MP*math.exp(-1*(1-Beta)*n*F*EeFe_Fe2p_MP/(R*TKelvin))));


Corrosionrate = iFe_Fe2p_MP * MWFe/ (n*F)

#Adjusting Ce_MP for the Ecorr
Ce_MP = math.exp((Ecorr + DeltaGFe_Fe2_MP/(n*F) + math.log(10)*2*R*TKelvin/(n*F)*pH)*(n*F)/(R*TKelvin))

#Solubility at M-P: CMPsat assuming that the solubility at the M_P interface is summation of what is diffused at this interface and the equilibrium concentration considering the developed Ecorr

CMP_Diffused = (0.476 * (1.101 + PhiOX) * ChiOX * DeltaOX) * Corrosionrate /(PhiOX * DFe * RhoOX * (1 - PhiOX))

CSat_MP = CMP_Diffused +  Ce_MP;


#Dissolution of magnetite at O_B interface Fe3O4 + 2H2O + 2H+ + 2e- = 3 Fe(OH)2 

EeFe2p_Fe3O4_OB = -1*DeltaGFe2p_Fe3O4_OB/(n*F) - math.log(10)*2*R*TKelvin/(n*F)*pH - 3*R*TKelvin/(n*F)*math.log(Ce_OB); # Equilibrium potential for dissolution reaction at O_B interface
i0Fe2p_Fe3O4_OB = F*kboltzman *TKelvin/h*(math.exp(-1*DeltaGFe2p_Fe3O4_OB/(R*TKelvin))*Ce_OB*math.exp(Beta*n*F*EeFe2p_Fe3O4_OB/(R*TKelvin)));
iFe3O4_Fe2p_OB = i0Fe2p_Fe3O4_OB *(math.exp(Beta*n*F/(R*TKelvin)*(Emixed-EeFe2p_Fe3O4_OB)) - math.exp(-1*(1-Beta)*n*F/(R*TKelvin)*(Emixed-EeFe2p_Fe3O4_OB)));

#Hydrogen consumption at O_B interface H2 = 2H+ + 2e- 
#Standard DeltaGH2_Hp is zero at any temperature

EeH2_Hp_OB =  - math.log(10)*2*R*TKelvin/(n*F)*pH - R*TKelvin/(n*F)*math.log(CeH2_OB); # At O_B interface
i0H2_Hp_OB = F*kboltzman *TKelvin/h*CeH2_OB*math.exp(Beta*n*F*EeH2_Hp_OB/(R*TKelvin));
iH2_Hp_OB = i0H2_Hp_OB *(math.exp(Beta*n*F/(R*TKelvin)*(Emixed-EeH2_Hp_OB)) - math.exp(-1*(1-Beta)*n*F/(R*TKelvin)*(Emixed-EeH2_Hp_OB)));

# Mixed Potential at the O_B interface
Emixed = R*Tkelvin/(Beta-(1-Beta)*n*F)*math.log((i0H2_Hp_OB*math.exp(Beta*n*F*EeH2_Hp_OB/(R*TKelvin))+i0Fe3O4_Fe2p_OB*math.exp(Beta*n*F*EeFe2p_Fe3O4_OB/(R*TKelvin)))/(i0H2_Hp_OB*math.exp(-1*(1-Beta)*n*F*EeH2_Hp_OB/(R*TKelvin))+i0Fe3O4_Fe2p_OB*math.exp(-1*(1-Beta)*n*F*EeFe3O4_Fe2p_OB/(R*TKelvin))));

#Adjusting the dissolution rate constant of magnetite
ked_OB = kd_OB *math.exp((1-Beta)*n*F*Emixed/(R*TKelvin))

#Adjusting the concentration of iron species at O_B interface 
Ce_OB = math.exp((Emixed + DeltaGFe3O4_Fe2_OB / (n * F) + 2*R*TKelvin*math.log(10)*pH/(n*F))*(-1)*n*F/(3*R*TKelvin))
CSat_OB = Ce_OB;

#H2 concentration

CH2generated = 0.005659 * Corrosionrate * (7.32 - PhiOX )


# H2 concentration at M_P interface assuminghigh mass transfer to the bulk

CeH2_MP = CH2coolant + (0.005659 * Corrosionrate * (7.32 - PhiOX )*(h*DeltaOX * ChiOX / (RhoOX * (1 - PhiOX)) + (DH2 * PhiOX))) / (DH2 * hH2 * PhiOX)

# H2 concentration at O_B interface

CeH2_OB = CH2coolant + (0.005659 * Corrosionrate * (7.32 - PhiOX ) / hH2)
 
 #Diffusivity of H2
 
DH2 = 0.0000222 * TKelvin *math.exp(-12400 / (R * T))


