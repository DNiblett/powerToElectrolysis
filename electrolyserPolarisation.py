#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 15:00:38 2025

@author: Daniel Niblett, Senior Research Associate, Newcastle University
"""

import numpy as np
import matplotlib.pyplot as plt
from types import SimpleNamespace

electrolyserTypes = ["PEM", "AEM", "Alkaline"]

plt.figure()  # one figure for all types

for electrolyserType in electrolyserTypes:
    print(f"Running simulation for: {electrolyserType}")

    # fresh parameter containers each loop
    p = SimpleNamespace(
        anode=SimpleNamespace(),
        cathode=SimpleNamespace(),
        electrolyte=SimpleNamespace()
    )
    cell = SimpleNamespace()

    # ---------------- Params per type ----------------
    if electrolyserType == "PEM":
        print("Using PEM parameters")
        p.anode.j0, p.cathode.j0 = 1, 1
        p.anode.b,  p.cathode.b  = 0.03, 0.04
        p.anode.CLthickness, p.cathode.CLthickness = 10e-6, 10e-6
        p.anode.specificArea, p.cathode.specificArea = 1e7, 1e7
        p.anode.E0, p.cathode.E0 = 1.23, 0
        p.electrolyte.thickness, p.electrolyte.thickness2 = 50e-6, 0
        p.electrolyte.conductivity, p.electrolyte.volFraction = 4, 1

  # Experimental data from https://doi.org/10.1016/j.jpowsour.2018.06.078
        experimentalData = [
            (0.0443223830146825, 1.377300095067075),
            (0.09124024807979603, 1.412992500264075),
            (0.17352306508322135, 1.4454843139326081),
            (0.25290256379302534, 1.4708144079433823),
            (0.3293546001901341, 1.4930812295341713),
            (0.4146433475682445, 1.5153163620999262),
            (0.49404095429235373, 1.5375726206823703),
            (0.5734566690307684, 1.5567550438364846),
            (0.6558119180914151, 1.5769515157916976),
            (0.7322880985075978, 1.5951198901447132),
            (0.8117038132460125, 1.6143023132988275),
            (0.8940650983114274, 1.633474173444597),
            (0.9675836363910727, 1.6537023344248443),
            (1.0528965277882572, 1.6718390197528257),
            (1.1323243145362083, 1.6889722192880532),
            (1.205848888620622, 1.7081757684588572),
            (1.2970468846670387, 1.7273159395795923),
            (1.3676621044530628, 1.7403823809020809),
            (1.4529508518311727, 1.7626175134678357),
            (1.5323726025743563, 1.7807753248125064),
            (1.6117943533175394, 1.7989331361571776),
            (1.6882826057432587, 1.8150522868913066),
            (1.7706559628182106, 1.8321749234181894),
            (1.8589204605471634, 1.8492764339283827),
            (1.9236143596553443, 1.8674870603147777),
            (0.03269100182589135, 1.351727051864371)
        ]

    elif electrolyserType == "AEM":
        print("Using AEM parameters")
        p.anode.j0, p.cathode.j0 = 0.0938, 1
        p.anode.b,  p.cathode.b  = 0.05, 0.04
        p.anode.CLthickness, p.cathode.CLthickness = 50e-6, 50e-6
        p.anode.specificArea, p.cathode.specificArea = 2e6, 2e6
        p.anode.E0, p.cathode.E0 = 1.23, 0
        p.electrolyte.thickness, p.electrolyte.thickness2 = 70e-6, 0
        p.electrolyte.conductivity, p.electrolyte.volFraction = 6, 1

# Exp data from: https://www.fuelcellstore.com/reinforced-aem-based-ccm
        experimentalData = [
            (0.022706852028724078, 1.5041152263374489),
            (0.04000141497753723, 1.5339506172839508),
            (0.07965262301460962, 1.565843621399177),
            (0.10032898227740639, 1.579218106995885),
            (0.20022639640595702, 1.6327160493827164),
            (0.300088436096077, 1.6759259259259263),
            (0.4016413739431886, 1.7109053497942388),
            (0.501475114082564, 1.7458847736625516),
            (0.6012840921150377, 1.7736625514403295),
            (0.7010930701475113, 1.8014403292181072),
            (0.8008985107361419, 1.8281893004115228),
            (0.9006968764370864, 1.8528806584362143),
            (1.0022073649580814, 1.8755144032921813),
            (1.1054405886306555, 1.8991769547325106),
            (1.2035126817361774, 1.9218106995884776)
        ]

    elif electrolyserType == "Alkaline":
        print("Using Alkaline parameters")
        p.anode.j0, p.cathode.j0 = 0.0938, 0.695
        p.anode.b,  p.cathode.b  = 0.05, 0.0828
        p.anode.CLthickness, p.cathode.CLthickness = 10e-6, 10e-6
        p.anode.specificArea, p.cathode.specificArea = 1e6, 1e6
        p.anode.E0, p.cathode.E0 = 1.23, 0
        p.electrolyte.thickness, p.electrolyte.thickness2 = 500e-6, 0
        p.electrolyte.volFraction = 0.5
        p.electrolyte.conductivity = 56 * p.electrolyte.volFraction**1.5

# Alkaline data from: https://doi.org/10.1016/j.electacta.2024.145161
        experimentalData = [
            (0.04074074074074073, 1.9472846252445521),
            (0.06111111111111113, 1.9891418741145515),
            (0.08271604938271605, 2.0200679124783556),
            (0.1037037037037037, 2.048262834783782),
            (0.15493827160493828, 2.1037250669005374),
            (0.20864197530864198, 2.148253839753536),
            (0.26049382716049385, 2.186410758056174),
            (0.31296296296296294, 2.2227450583552586),
            (0.36543209876543215, 2.2545256245924126),
            (0.41728395061728396, 2.2835750747711887)
        ]

    else:
        print("Unknown electrolyser type; skipping.")
        continue

    cell.area = 1.0
    Asa = p.anode.CLthickness * p.anode.specificArea
    Asc = p.cathode.CLthickness * p.cathode.specificArea

    # model residual and derivative
    f  = lambda b_a, j, j0_a, L, k, b_c, j0_c, phi_anode, E0a, Asa, Asc: (
         b_a*np.log(j/(j0_a*Asa)) + (j*L)/k + b_c*np.log(j/(j0_c*Asc)) - phi_anode + E0a)
    dfdj = lambda b_a, j, L, k, b_c: (b_a/j) + (L/k) + (b_c/j)

    # sweep potentials
    dx = 0.01
    potential_values = np.arange(1.23, 3.0, dx)

    currentDensityResult = []
    potentialResult = []

    for phi_anode in potential_values:
        j = 0.5  # initial guess (A/m^2)
        for _ in range(20):
            f1  = f(p.anode.b, j, p.anode.j0, p.electrolyte.thickness, p.electrolyte.conductivity,
                    p.cathode.b, p.cathode.j0, phi_anode, p.anode.E0, Asa, Asc)
            df1 = dfdj(p.anode.b, j, p.electrolyte.thickness, p.electrolyte.conductivity, p.cathode.b)
            j  -= f1/df1
            if j <= 0:  # guard against negative/zero
                j = 1e-6
        currentDensityResult.append(j)
        potentialResult.append(phi_anode)

    currentDensityResult = np.array(currentDensityResult)
    potentialResult     = np.array(potentialResult)

    color = {'PEM': 'k', 'AEM': 'g', 'Alkaline': 'b'}.get(electrolyserType, 'r')

    # model curve
    plt.plot(currentDensityResult/1e4, potentialResult, color=color, linewidth=0.7,
             label=f"Model {electrolyserType}")

    # experimental markers
    exp = np.array(experimentalData)
    plt.plot(exp[:, 0], exp[:, 1], marker='o', linestyle='none',
             markerfacecolor=color, markeredgecolor='k',
             label=f"Exp {electrolyserType}")

# --------------- Figure formatting (once) ---------------
plt.xlabel("Current Density (A/cm²)")
plt.ylabel("Potential (V)")
plt.title("Electrolyser Polarisation Curve")
plt.grid(True)
plt.xlim(0, 2)
plt.ylim(1.23, 3)
plt.legend()
plt.tight_layout()
plt.show()
