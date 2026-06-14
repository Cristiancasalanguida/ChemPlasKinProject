"""
Code for solving mass flow rate in a CRN system.
"""

import numpy as np
import scipy as sp
from pathlib import Path

# SICCA FUNCTIONS
def solveFull(air_in, fuel_in, parameters, ox_Y, fuel_Y, f_st):
    """
    Solves for the mass flow rate in a CRN system given input parameters.

    Parameters:
    air_in (float): Mass flow rate of air input.
    fuel_in (float): Mass flow rate of fuel input.
    parameters (list): List containing necessary parameters for the system.
    ox_Y (dict): Oxidizer composition dictionary for every reactor. (e.g., {'PSR1': 0.12, 'PSR2': 0.4})
    fuel_Y (dict): Fuel composition dictionary for every reactor. (e.g., {'PSR1': 0.88, 'PSR2': 0.01})
    f_st (float): Stoichiometric fuel-to-air ratio.

    Returns:
    massFlowRate (numpy.ndarray): The calculated mass flow rate vector.
    """
    p1, p2, p3, p4, p5, p6, p7, p8, p9, p10 = parameters

    mair_I      = air_in * p1 * (1 - p3) / (1 + p1)
    mair_II     = air_in * (1 - p3) / (1 + p1)
    mair_III    = air_in * p3
    mfuel_I     = fuel_in * p2 / (1 + p2)
    mfuel_II    = fuel_in / (1 + p2)

    mOut = air_in + fuel_in

    yO2_air_II = ox_Y['airInlet']
    yH2_r1 = fuel_Y['PSR1']
    yO2_r1 = ox_Y['PSR1']
    yH2_r2 = fuel_Y['PSR2']
    yO2_r2 = ox_Y['PSR2']
    yH2_r3 = fuel_Y['PSR3']
    yO2_r3 = ox_Y['PSR3']
    yH2_r4 = fuel_Y['PSR4']
    yO2_r4 = ox_Y['PSR4']

    cA =  (yH2_r1 - f_st * yO2_r1)
    cB = -(yH2_r2 - f_st * yO2_r2)
    cC = -(yH2_r2 - f_st * yO2_r2)
    cD =  (yH2_r3 - f_st * yO2_r3)
    cE = -(yH2_r2 - f_st * yO2_r2)
    cF =  (yH2_r4 - f_st * yO2_r4)

    matA = np.array([ # (13, 13) matrix for the 13 unknown mass flow rates
        [ 1, -1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # mass conservation r1 (out=positive)
        [cA, cB, cC, cD, cE, cF,  0,  0,  0,  0,  0,  0,  0], # mass conservation r2 with phi=1
        [ 0,  0,  1, -1,  0,  0, -1,  1, -1,  1,  0,  0,  0], # mass conservation r3 (in=positive)
        [ 0,  0,  0,  0,  1, -1,  0,  0,  1, -1, -1,  1,  0], # mass conservation r4 (in=positive)
        [ 0,  0,  0,  0,  0,  0, -1, 1,  0,   0, -1,  1,  1], # mass conservation r5 (out=positive)
        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1], # mass conservation r6 (in=positive)
        [ 0,  0, -1,  0, p4,  0,  0,  0,  0,  0,  0,  0,  0], # mC = p4*mE
        [ 0,  0,  0, p5,  0,  0, -1,  0,  0,  0,  0,  0,  0], # mG = p5*mD
        [ 1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0], # mA = p6*mair_II
        [ 0,  0,  0,  0,  0,  0,  0, p7,  0,  0,  0,  0, -1], # mO = p7*mH
        [ 0,  0,  0,  0,  0, -1,  0,  0,  0, p8,  0,  0,  0], # mF = p8*mL
        [ 0,  0,  0,  0,  0,  0, -1,  0, p9,  0,  0,  0,  0], # mG = p9*mI
        [ 0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0]  # mN = p10*mair_III
    ])

    
    vecb = np.array([ # (13,) vector for the 13 unknown mass flow rates
        (mair_I  + mfuel_I), # mass flow inlet in r1
        (f_st*yO2_air_II*mair_II - mfuel_II), # mass flow inlet in r2 with phi=1
        0, # mass flow inlet in r3
        0, # mass flow inlet in r4
        mair_III, # mass flow inlet in r5
        air_in + fuel_in, # mass flow outlet in r6
        0, # p4 constraint
        0, # p5 constraint
        p6*mair_II, # p6 constraint
        0, # p7 constraint
        0, # p8 constraint
        0, # p9 constraint
        p10*mair_III # p10 constraint
    ])
    
    m_alphabet = sp.linalg.solve(matA, vecb) # (,13) vector of the 13 unknown mass flow rates
    m_alphabet = m_alphabet.reshape(-1,) # ensure it's a 1D array for easier unpacking
    flow_rates = np.concatenate((
        np.array([mair_I, mair_II, mair_III, mfuel_I, mfuel_II, mOut]),
        m_alphabet,
    )) # (19,) vector: 6 pre-calculated + 13 unknowns from m_alphabet

    # Column order of m_alphabet: mA(0) mB(1) mC(2) mD(3) mE(4) mF(5) mG(6) mH(7) mI(8) mL(9) mM(10) mN(11) mO(12)
    names = ['airI', 'airII', 'airIII', 'fuelI', 'fuelII', 'mOut',
             'mA', 'mB', 'mC', 'mD', 'mE', 'mF', 'mG', 'mH', 'mI', 'mL', 'mM', 'mN', 'mO']
    
    massFlows = {}
    for i, mfr in enumerate(flow_rates):
        massFlows[names[i]] = mfr
    return massFlows

def mainSolver2(ox_Y, fuel_Y):
    air_in = 0.5
    fuel_in = 0.05
    f_st = 1/34.3
    parameters = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1] # Placeholder parameters, replace with actual values
    massFlowRates = solveFull(air_in, fuel_in, parameters, ox_Y, fuel_Y, f_st)
    return massFlowRates

#FOLCARELLI ET AL. FUNCTIONS
def mainSolverFOLC(ox_Y, fuel_Y):
    f_st_air = 1/34.3
    f_st_O2 = 1/8.0
    fuel_in = 0.048e-3 # kg/s
    phi_global = 0.17
    f = phi_global * f_st_air
    air_in = fuel_in / f
    p1, p2, p3, p4 = [0.189, 0.125, 0.216, 0.099] # Naik
    # p1, p2, p3, p4 = [0.172, 0.125, 0.216, 0.099] # Creck
    # p1, p2, p3, p4 = [0.217, 0.125, 0.216, 0.099] # Gotama

    yH2_r1 = fuel_Y['PSR1']
    yO2_r1 = ox_Y['PSR1']
    yH2_r2 = fuel_Y['PSR2']
    yO2_r2 = ox_Y['PSR2']
    yH2_r3 = fuel_Y['PSR3']
    yO2_r3 = ox_Y['PSR3']
    yH2_r5 = fuel_Y['PSR5']
    yO2_r5 = ox_Y['PSR5']

    mAir_I = air_in * p1 / (1 + p1)
    mAir_II = air_in - mAir_I
    mH2 = fuel_in

    mA = fuel_in / (1 + p3)
    mB = p2 * mAir_I #!!!!!!!!!!!!!!!!!
    mC = mAir_I * (1 + p2)
    mF = p3 * mA
    mI = air_in + fuel_in
    mG = p4 * mI
    mH = mG + mF # + mI - mA - mC +mB - mAir_II

    # !!!!!!!!! NEL PAPER MANCA IL TERMINE mC * yH2_r3
    K1 = mA * yH2_r1 + mC * yH2_r3 - (mA + mC) * yH2_r2 # mA * yH2_r1 - (mB + mD) * yH2_r2
    # !!!!!!!!!
    K2 = mA * yO2_r1 + mC * yO2_r3 - (mA + mC) * yO2_r2 # mA * yO2_r1 + mC * yO2_r3 - (mB + mD) * yO2_r2
    K3 = (f_st_O2 * yO2_r5 - yH2_r5) - (f_st_O2 * yO2_r2 - yH2_r2)
    try:
        mE = (K1 - f_st_O2 * K2) / K3
    except ZeroDivisionError:
        mE = 0.0

    mE = max([mE, 0.0])

    mD = mA + mC + mE -mB

    mOut = mI

    massFlows = {
        'airI': mAir_I,
        'airII': mAir_II,
        'fuelI': mH2,
        'out': mOut,
        'mA': mA,
        'mB': mB,
        'mC': mC,
        'mD': mD,
        'mE': mE,
        'mF': mF,
        'mG': mG,
        'mH': mH,
        'mI': mI
    }

    return massFlows

#MADONIA ET AL. FUNCTIONS
def mainSolver(ox_Y, fuel_Y):
    params_file = Path(__file__).parent / "massflow_params.txt"
    with open(params_file) as f:
        lines = f.readlines()
    air_in     = float(lines[0].strip())
    phi_global = float(lines[1].strip())

    f_st_air = 1/34.3
    f_st_O2 = 1/8.0

    parameters = [0.57, 0.306, 0.708, 0.293, 0.141, 0.960, 0.278, 0.142, 0.130] # theta1 to theta9 Naik
    # parameters = [0.572, 0.312, 0.718, 0.280, 0.146, 0.912, 0.288, 0.148, 0.129] # theta1 to theta9 Creck
    
    theta1, theta2, theta3, theta4, theta5, theta6, theta7, theta8, *_ = parameters

    air_II = air_in / (theta2 + 1) # kg/s
    air_I = theta2 * air_II  # kg/s
    phi2 = phi_global * air_in / (theta1 * air_I + air_II)
    phi1 = theta1 * phi2
    m1 = (1 + f_st_air * phi1) * air_I  # kg/s
    m2 = (1 + f_st_air * phi2) * air_II  # kg/s
    m2ORZ = theta3 * m2  # kg/s
    m2FLSW = m2 - m2ORZ  # kg/s

    yH2r1 = (phi1 * f_st_air) / (1 + phi1 * f_st_air)
    yO2r1 = (1 - yH2r1) * 0.219 * 32 / 29

    yH2r2 = (phi2 * f_st_air) / (1 + phi2 * f_st_air)
    yO2r2 = (1 - yH2r2) * 0.219 * 32 / 29

    zH2_ZAX = fuel_Y['ZAX']
    zO2_ZAX = ox_Y['ZAX']
    zH2_ORZ = fuel_Y['ORZ']
    zO2_ORZ = ox_Y['ORZ']
    zH2_FLSW = fuel_Y['FLSW']
    zO2_FLSW = ox_Y['FLSW']
    zH2_FLAX = fuel_Y['FLAX']
    zO2_FLAX = ox_Y['FLAX']
    zH2_IRZ = fuel_Y['IRZ']
    zO2_IRZ = ox_Y['IRZ']
    zH2_res2 = yH2r2
    zO2_res2 = yO2r2

    denThetas = (theta5 * theta6 * theta7 * theta8 
                - theta5 * theta6 * theta8
                - theta5 * theta7 * theta8
                + theta5
                - theta6 * theta7 * theta8
                + theta6 * theta8
                + theta7 * theta8
                -1)
    
    zres2 = zH2_res2 - f_st_O2 * zO2_res2 * phi_global
    zZAX = zH2_ZAX - f_st_O2 * zO2_ZAX * phi_global
    zIRZ = zH2_IRZ - f_st_O2 * zO2_IRZ * phi_global
    zORZ = zH2_ORZ - f_st_O2 * zO2_ORZ * phi_global
    zFLSW = zH2_FLSW - f_st_O2 * zO2_FLSW * phi_global

    z1 = zres2 - zFLSW
    z2 = zZAX - zFLSW
    z3 = zIRZ - zFLSW
    z4 = zORZ - zFLSW
    kE = (1 - theta5) * theta6
    kB = theta4 / (1 - theta4)
    kA = (1 - theta8) / denThetas
    kT = (m1 + theta8 * m2) / denThetas

    numA = z3 * kE * kT + z4 * kB * kE * kT - z1 * m2FLSW - z4 * (1 + kB) * m2ORZ - z4 * kB * m2FLSW
    denA = z2 + z3 * kE * kA + z4 * kB + z4 * kB * kE * kA

    mA = (numA / denA) if abs(denA) > 1e-20 else 0.0
    mA_max1 = (m1 - theta5 * kT) / (1 - theta5 * kA)  # ZAX balance: guarantees mF >= 0 and TFH >= 0
    mA_max2 = (m1 + theta8 * m2) / (1 - theta8)
    mA_max = min(mA_max1, mA_max2)
    mA = min(max(0.0, mA), mA_max)
    TFH = kA * mA - kT  # = mF + mH >= 0 by construction

    mE = kE * TFH
    mD = m2 + mA + mE
    mB = kB * mD
    mC = mB + m2ORZ
    mI = (1 - theta5) * (1 - theta6) * theta7 * TFH
    mH = theta8 * (mD + mI)
    mG = theta5 * TFH
    mF = m1 + mG - mA
    mM = (1 - theta5) * (1 - theta6) * (1 - theta7) * TFH  # Eq. 23i — IRZ balance
    mN = (1 - theta8) * (mD + mI)                           # Eq. 23j
    mOut = mM + mN

    # massflows = np.array([m1, m2FLSW, m2ORZ, mA, mB, mC, mD, mE, mF, mG, mH, mI, mM, mN, mOut])

    massFlows = {
        'm1': m1,
        'm2FLSW': m2FLSW,
        'm2ORZ': m2ORZ,
        'mA': mA,
        'mB': mB,
        'mC': mC,
        'mD': mD,
        'mE': mE,
        'mF': mF,
        'mG': mG,
        'mH': mH,
        'mI': mI,
        'mM': mM,
        'mN': mN
    }

    return massFlows

#TEST FUNCTION
def mainSolverSERIAL(ox_Y, fuel_Y):
    f_st = 1/34.3
    fuel_in = 0.005e-3 # kg/s
    phi_global = 1.5
    f = phi_global * f_st
    air_in = fuel_in / f
    
    mAir_I = air_in
    mH2 = fuel_in

    mA = fuel_in + air_in

    mOut = mA
    
    massFlows = {
        'airI': mAir_I,
        'fuelI': mH2,
        'out': mOut,
        'mA': mA,
    }
    """
    massFlows = {
        'airI': 0.0,
        'fuelI': 0.0,
        'out': 0.0,
        'mA': 0.0,
    }
    """

    return massFlows
