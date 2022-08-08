import numpy as np

# define useful values
R = 8.31434  # J/mol.K
F = 96484.56  # C/mol

pO2 = 1  # partial pressures, both in atm
pH2 = 1  # for now =1

#at 298k
delta_G_nio3 = -541800
delta_G_nip2 = -46400
delta_G_nio2 = -453100
delta_G_w = -237178
delta_G_nin6p2 = -250000
delta_G_nin4p2 = -192000
delta_G_o2a = 16300
delta_G_n = -26600
delta_G_ng = -16500
delta_G_nhp = -79370


def R1(x, T): #NI(2+)+ 2e- = Ni
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltah1 = 64000
    delta_G_R1_0 = -delta_G_nip2
    delta_G_R1 = (T*delta_G_R1_0/298)+T*deltah1*((1/T)-(1/298))
    E_theta_1 = -delta_G_R1/(2*F)
    V1 = E_theta_1 + (C2/2)*np.log10(x)
    return V1


def R2(pH, T): #Ni(OH)2 +2H+ +2e- =Ni + 2H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltah2 = -33560
    delta_G_R2_0 = 2*delta_G_w - delta_G_nio2
    delta_G_R2 = (T*delta_G_R2_0/298)+T*deltah2*((1/T)-(1/298))
    E_theta_2 = -delta_G_R2/(2*F)
    V2 = []
    if type(pH) == list:  # for list(vector), below is only for single value
        for pHval in pH:
            V2.append(E_theta_2 - (C2*2*pHval)/2)
    else:
        V2 = E_theta_2 - (C2*2*pH)/2
    return V2


def R3(pH, nhp, nin4p2, T): #Ni(NH3)4(2+) + 4H+ +2e- = Ni + 4NH4+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltah3 = -530000
    delta_G_R3_0 = 4*delta_G_nhp - delta_G_nin4p2
    delta_G_R3 = (T*delta_G_R3_0/298)+T*deltah3*((1/T)-(1/298))
    E_theta_3 = - delta_G_R3/(2*F)
    V3 = []
    for pHval in pH:
        V3.append(E_theta_3 - (C2*4*pHval/2) - (C2*4*np.log10(nhp)/2) + (C2*np.log10(nin4p2)/2))
    return V3


def R4(pH, nhp, nin6p2, T): #Ni(NH3)+6H+ +2e-=Ni +6NH4+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltah4 = -795000
    delta_G_R4_0 = 6*delta_G_nhp - delta_G_nin6p2
    delta_G_R4 = (T*delta_G_R4_0/298)+T*deltah4*((1/T)-(1/298))
    E_theta_4 = -delta_G_R4/(2*F)
    V4 = []
    for pHval in pH:
        V4.append(E_theta_4 - (C2*6*pHval/2) - (C2*6*np.log10(nhp)/2) + (C2*np.log10(nin6p2)/2))
    return V4


def R5(n, nin4p2, T): #Ni(NH3)4(2+) +2e- = Ni +4NH3
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltah5 = -321160
    delta_G_R5_0 = 4*delta_G_n - delta_G_nin4p2
    delta_G_R5 = (T*delta_G_R5_0/298)+T*deltah5*((1/T)-(1/298))
    E_theta_5 = -delta_G_R5/(2*F)
    V5 = E_theta_5 - (C2*4*np.log10(n)/2) + (C2*np.log10(nin4p2)/2)
    return V5


def R6(n, nin6p2, T): # Ni(NH3)+6H+ +2e- =Ni + 6NH3
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltah6 = -481740
    delta_G_R6_0 = 6*delta_G_n - delta_G_nin6p2
    delta_G_R6 = (T*delta_G_R6_0/298)+T*deltah6*((1/T)-(1/298))
    E_theta_6 = -delta_G_R6/(F*2)
    V6 = E_theta_6 - (C2*6*np.log10(n)/2) + (C2*np.log10(nin6p2)/2)
    return V6


def V1(x, T): #Ni(OH)2 +2H= Ni+ 2H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltahv1 = -33560
    delta_G_V1_0 = 2*delta_G_w + delta_G_nip2 - delta_G_nio2
    delta_G_V1 = (T*delta_G_V1_0/298)+T*deltahv1*((1/T)-(1/298))
    logK1 = delta_G_V1/ C1
    pH_V1 = (logK1 - np.log10(x))/2
    return pH_V1


def V2(nhp, nin4p2, T): #Ni(OH)2 +4NH4+ = Ni(NH3)4(2+) +2H2O +2H+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltahv2 = 496440
    delta_G_V2_0 = 2*delta_G_w + delta_G_nin4p2 - delta_G_nio2 - 4*delta_G_nhp
    delta_G_V2 = (T*delta_G_V2_0/298)+T*deltahv2*((1/T)-(1/298))
    logK2 = delta_G_V2/C1
    pH_V2 = (18.26 + np.log10(nin4p2) - 4*np.log10(nhp))/2
    return pH_V2


def V3(nhp, nin4p2, nin6p2,T): # Ni(NH3)4(2+) + 2NH4+ = Ni(NH3)6(2+) +2H+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltahv3 = 265000
    delta_G_V3_0 = delta_G_nin6p2 - 2*delta_G_nhp - delta_G_nin4p2
    delta_G_V3 = (T*delta_G_V3_0/298)+T*deltahv3*((1/T)-(1/298))
    logK3 = delta_G_V3/C1
    pH_V3 = (np.log10(nin6p2) - np.log10(nin4p2) - logK3 - 2*np.log10(nhp))/2
    return pH_V3


def V4(): #NH4+ =H+ +NH3
    delta_G_V4_0 = 52806  # calculated using 1st law
    logK4 = 0
    pH_V4 = 9.25
    return pH_V4


def V5(n, nin6p2, T): #Ni(OH)2 +6NH3 +2H+ = Ni(Nh3)4(2+) +2H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltahv5 = 448180
    delta_G_V5_0 = 2*delta_G_w + delta_G_nin6p2 - delta_G_nio2 - 6*delta_G_n
    delta_G_V5 = (T*delta_G_V5_0/298)+T*deltahv5*((1/T)-(1/298))
    logK5 = delta_G_V5/C1
    pH_V5 = (logK5 + 6*np.log10(n) - np.log10(nin6p2))/2
    return pH_V5


def T1(nip2, pH, T): #Ni(OH)3 +3H+ +e- =Ni(2+) +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltaht1 = -243290
    delta_G_T1_0 = 3*delta_G_w + delta_G_nip2 - delta_G_nio3
    delta_G_T1 = (T*delta_G_T1_0/298)+T*deltaht1*((1/T)-(1/298))
    E_theta_T1 = -delta_G_T1/(F)
    if type(pH) == list:
        vt1 = []
        for pHval in pH:
            vt1.append(E_theta_T1 - C2*np.log10(nip2) - C2*3*pHval)
    else:
        vt1 = E_theta_T1 - C2*np.log10(nip2) - C2*3*pH
    return vt1


def T2(pH, T): #Ni(OH)3 +H+ +e- =Ni(OH)2 +H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltaht2 = -145730
    delta_G_T2_0 = delta_G_nio2 + delta_G_w - delta_G_nio3
    delta_G_T2 = (T*delta_G_T2_0/298)+T*deltaht2*((1/T)-(1/298))
    E_theta_T2 = -delta_G_T2/F
    if type(pH) == list:
        vt2 = []
        for pHval in pH:
            vt2.append(E_theta_T2 - (C2*1*pHval)/1)
    else:
        vt2 = E_theta_T2 - (C2*1*pH)/1
    return vt2


def T3(nhp, nin4p2, pH, T): #Ni(OH)3 +4NH4+ +e- =Ni(Nh3)4(2+) +3H2O +H+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltaht3 = 350710
    delta_G_T3_0 = 3*delta_G_w + delta_G_nin4p2 - 4*delta_G_nhp - delta_G_nio3
    delta_G_T3 = (T*delta_G_T3_0/298)+T*deltaht3*((1/T)-(1/298))
    E_theta_T3 = - delta_G_T3/F
    if type(pH) == list:
        vt3 = []
        for pHval in pH:
            vt3.append(E_theta_T3 + C2*pHval + C2*4*np.log10(nhp) - C2*np.log10(nin4p2))
    else:
        vt3 = E_theta_T3 + C2*pH + C2*4*np.log10(nhp) - C2*np.log10(nin4p2)
    return vt3


def T4(pH, nhp, nin6p2, T): #Ni(OH)3 +6NH4+ +e- =Ni(Nh3)6(2+) +3H2O +H+
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltaht4 = 615710
    delta_G_T4_0 = 3*delta_G_w + delta_G_nin6p2 - 6*delta_G_nhp - delta_G_nio3
    delta_G_T4 = (T*delta_G_T4_0/298)+T*deltaht4*((1/T)-(1/298))
    E_theta_T4 = -delta_G_T4/F
    if type(pH) == list:
        vt4 = []
        for pHval in pH:
            vt4.append(E_theta_T4 + C2*3*pHval + C2*6*np.log10(nhp) - C2*np.log10(nin6p2))
    else:
        vt4 = E_theta_T4 + C2*3*pH + C2*6*np.log10(nhp) - C2*np.log10(nin6p2)
    return vt4


def T5(n, pH, nin6p2, T): #Ni(OH)3 +6NH3 +3H +e- =Ni(Nh3)6(2+) +3H2O
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltaht5 = 302450
    delta_G_T5_0 = 3*delta_G_w + delta_G_nin6p2 - 6*delta_G_n - delta_G_nio3
    delta_G_T5 = (T*delta_G_T5_0/298)+T*deltaht5*((1/T)-(1/298))
    E_theta_T5 = -delta_G_T5/F
    if type(pH) == list:
        vt5 = []
        for pHval in pH:
            vt5.append(E_theta_T5 + 6*C2*np.log10(n) - 3*C2*pHval - C2*np.log10(nin6p2))
    else:
        vt5 = E_theta_T5 + 6*C2*np.log10(n) - 3*C2*pH - C2*np.log10(nin6p2)
    return vt5

def W1(pH, T): #O2+ 4H+ +4e- =2H20
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    deltahw1 = -571660
    delta_G_W1_0 = 2*delta_G_w
    delta_G_W1 = (T*delta_G_W1_0/298)+T*deltahw1*((1/T)-(1/298))
    E_theta_W1 = -delta_G_W1/(4*F)
    VW1 = []
    for pHval in pH:
        VW1.append(E_theta_W1 - (C2*4*pHval/4) + (C2*1*np.log10(pO2)/4))
    return VW1


def W2(pH, T):   #2H+ 2e- = H2
    C1 = -np.log(10)*R*T
    C2 = -C1/F
    VW2 = []
    for pHval in pH:
        VW2.append(-(C2*2*pHval/2) - (C2*1/2))
    return VW2
