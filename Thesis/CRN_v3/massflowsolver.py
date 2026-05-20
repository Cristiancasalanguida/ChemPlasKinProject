"""
Code for solving mass flow rate in a CRN system.
"""

import numpy as np
import scipy as sp

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

def mainSolver(ox_Y, fuel_Y):
    air_in = 0.5
    fuel_in = 0.05
    f_st = 1/8
    parameters = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1] # Placeholder parameters, replace with actual values
    massFlowRates = solveFull(air_in, fuel_in, parameters, ox_Y, fuel_Y, f_st)
    return massFlowRates