#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  6 16:26:30 2025

Electrolyser Model

@author: Dr Daniel Niblett - Newcastle University
"""

"""
Electrolyser Model for PEM, Alkaline, AEM and Membraneless:

2. Kinetic overpotential from rate of reaction at cathode (HER)
3. Ionic conduction through electrolyte (solid, fluid or porous)
5. Kinetic overpotential from rate of reaction at anode (OER)
6. Cell potential is potential of anode BPP - cathode BPP.

Assumptions: Isothermal, 0D, steady-state, single-phase, thin-catalyst layer
with surface roughness. No losses through catalyst layer.

"""


# Make this a function


# Example usage: runElectrolyser(powerInput=10, stackNumber=3, stackCapacity=2, minLoad=10, electrolyserType="PEM", runMode="transient")

def runElectrolyserFull(powerInput, stackNumber, stackCapacity, minLoad, electrolyserType, runMode):

    import numpy as np
    import matplotlib.pyplot as plt
    from types import SimpleNamespace
    
    p = SimpleNamespace();
    p.anode = SimpleNamespace();
    p.cathode = SimpleNamespace();
    p.electrolyte = SimpleNamespace();
    cell = SimpleNamespace();
    
    #mode = "validate-polarisation"
    mode = runMode;
    
    
    #powerInput = np.array([0,1,3,4,5])
    #stackNumber = 5;
    #stackCapacity = 1;
    #runMode = "transient"
    
    
    if mode == 'transient':
        
        
        numberOfStacks = stackNumber;
        stackCapacity = stackCapacity*1e6;
        
        
        powerInput = powerInput*1e6/numberOfStacks; # assume 3MW
       
        
        # new function to solve for j_cell
        # f = N * A_cell * V_cell(j_cell) * j_cell - power;
        # df = df/dj
        
        #electrolyserType = "PEM";
        
        """
         power = N * A * V * j
         set point of 1.8 V (j = 2 PEM)
         set point of 2.2 V (j = 0.25 Alkaline)
         2 MW system: 2e6/(2.2*2500) = N * A_cells
         A_cells = 0.5 * 0.5 [commercial MW stacks areas of 0.1 - 0.3 m2 https://doi.org/10.1021/acs.chemrev.3c00904]
         N_cells = 2e6/(2.2*2500*1x1) = N
         
        """
        
        

        
        # solve for current and output current
        
        # Read electrolyser type
        
        if electrolyserType == "PEM":
            # Do PEM setup
           # print("Using PEM parameters")
            
            # kinetic properties
            p.anode.j0 = 1;         # [A/m2]
            p.cathode.j0 = 1;       # [A/m2]
            p.anode.b = 0.03;       # [V]
            p.cathode.b = 0.04;     # [V]
            
            # electrode properties 
            p.anode.CLthickness = 10e-6; # [m]
            p.anode.specificArea = 1e+7; # [m2/m3]
            p.anode.E0 = 1.23; # [V]
            
            p.cathode.CLthickness = 10e-6; # [m]
            p.cathode.specificArea = 1e+7; # [m2/m3]
            p.cathode.E0 = 0; # [V]
            
            p.electrolyte.thickness = 50e-6; # [m]
            p.electrolyte.thickness2 = 0;    # [m] - additional layer thickness
            p.electrolyte.conductivity = 4; # [S/m]
            p.electrolyte.volFraction = 1;   # [S/m]
            
            
            cellArea = 0.5*0.5; # width x height
            numberOfCells = round(stackCapacity/(1.8*20000*cellArea)); # number of cells in stack
            
        
        elif electrolyserType == "AEM":
            # Do AEM setup
            #print("Using AEM parameters")
            # kinetic properties
            p.anode.j0 = 0.0938;         # [A/m2]
            p.cathode.j0 = 1;       # [A/m2]
            p.anode.b = 0.05;       # [V]
            p.cathode.b = 0.04;     # [V]
            
            # electrode properties 
            p.anode.CLthickness = 50e-6; # [m]
            p.anode.specificArea = 2e+6; # [m2/m3]
            p.anode.E0 = 1.23; # [V]
            
            p.cathode.CLthickness = 50e-6; # [m]
            p.cathode.specificArea = 2e+6; # [m2/m3]
            p.cathode.E0 = 0; # [V]
            
            p.electrolyte.thickness = 70e-6; # [m]
            p.electrolyte.thickness2 = 0;    # [m] - additional layer thickness
            p.electrolyte.conductivity = 6; # [S/m]
            p.electrolyte.volFraction = 1;   # [S/m]
            
            cellArea = 0.5*0.5; # width x height
            numberOfCells = round(stackCapacity/(1.8*6000*cellArea)); # number of cells in stack
        
            
        elif electrolyserType == "Alkaline":
            # Do AEM setup
            #print("Using AEM parameters")
            # kinetic properties
            p.anode.j0 = 0.0938;         # [A/m2]
            p.cathode.j0 = 0.695;       # [A/m2]
            p.anode.b = 0.05;       # [V]
            p.cathode.b = 0.0828;     # [V]
            
            # electrode properties 
            p.anode.CLthickness = 10e-6; # [m]
            p.anode.specificArea = 1e+6; # [m2/m3]
            p.anode.E0 = 1.23; # [V]
            
            p.cathode.CLthickness = 10e-6; # [m]
            p.cathode.specificArea = 1e+6; # [m2/m3]
            p.cathode.E0 = 0; # [V]
            
            p.electrolyte.thickness = 500e-6; # [m]
            p.electrolyte.thickness2 = 0;    # [m] - additional layer thickness
            p.electrolyte.volFraction = 0.5;   # [S/m]
            p.electrolyte.conductivity = 56*p.electrolyte.volFraction**(1.5); # [S/m]
            
            cellArea = 1.3*1.3; # width x height [working backwards from 230 cells https://nelhydrogen.com/wp-content/uploads/2024/08/A-Series-Spec-Sheet%E2%80%93DOC001974_03.pdf]
            numberOfCells = round(stackCapacity/(2.2*2500*cellArea)); # number of cells in stack
            
        else:
            print("Unknown electrolyser type")
        
        
        # Then 
        #update parameters
          
        # Cell Properties
        Asa = p.anode.CLthickness*p.anode.specificArea;
        Asc = p.cathode.CLthickness*p.cathode.specificArea;
        
        # Boundary Conditions
        phi_cathode = 0;
        #phi_anode = specifiedPotential;
        j = 100;
        result = [];
        potentialResult = [];
        
        #f = lambda b_a,j,j0_a,L,k,b_c,j0_c,phi_anode,E0a,Asa,Asc : b_a * np.log(j/(j0_a*Asa)) + j*L/k + b_c*np.log(j/(j0_c*Asc)) - phi_anode + E0a;
        #df = lambda b_a,j,L,k,b_c : b_a/j + L/k + b_c/j;
        
        f = (lambda b_a,j,j0_a,L,k,b_c,j0_c,E0a,Asa,Asc,A_cell,N_cells,power : 
             (b_a * np.log(j/(j0_a*Asa)) + j*L/k + b_c*np.log(j/(j0_c*Asc)) + E0a)*j*A_cell*N_cells - power);
            
        df = lambda b_a,j,j0_a,L,k,b_c,j0_c,E0a,Asa,Asc,A_cell,N_cells,power  : A_cell*N_cells * (
        b_a*np.log(j/(j0_a*Asa)) +
        b_c*np.log(j/(j0_c*Asc)) +
        E0a +
        (b_a + b_c) +
        2*j*L/k);
    #%%    
        power = powerInput;
        A_cell = cellArea;
        N_cells = numberOfCells;
        
        #potential = phi_anode;
        
        # Set up parameters
        iterations = 10
    
        
        #j = 0.5  # initial guess for j
        j = np.ones(len(powerInput))*1000;
        result = []
    
        for i in range(iterations):
            #phi_anode = potential  # Constant in this case
    
            f1 = f(
                p.anode.b, j, p.anode.j0,
                p.electrolyte.thickness, p.electrolyte.conductivity,
                p.cathode.b, p.cathode.j0, p.anode.E0, Asa, Asc,
                A_cell, N_cells, power
            )
    
            df1 = df(
                p.anode.b, j, p.anode.j0,
                p.electrolyte.thickness, p.electrolyte.conductivity,
                p.cathode.b, p.cathode.j0, p.anode.E0, Asa, Asc,
                A_cell, N_cells, power
            )
    
            jn = j - f1 / df1
            j = jn
            result.append(j)
        
  #%%          
        predictedCurrentDensity = result[i];
        
        
    print(j)
    

    # then get back the potential by the power equation P = I*V
    potential = power/(j*A_cell*N_cells);
    
    cellCurrent = predictedCurrentDensity*cellArea;
    hydrogenProductionRate_moles = (cellCurrent/(2*96485))*numberOfCells*numberOfStacks; # moles per second
    hydrogenProductionRate_kg = hydrogenProductionRate_moles*0.002; # kg/s
    H2_per_hour = hydrogenProductionRate_kg*3600;
    
    predictedEfficiency = 1.48/potential;
    powerUsed = potential*numberOfCells*cellCurrent*numberOfStacks/1e6; #MW
    powerUsedPerStack = powerUsed/numberOfStacks;
    
    #stackRatedPower = 20000*1.8*numberOfCells*cellArea/1e6; # stack rated power (MW)
    stackRatedPower = stackCapacity*(1e-6); # convert back to MW
    capacityPerStack = (powerUsedPerStack/stackRatedPower)*100;
    
    """
    if (capacityPerStack < minLoad):
        potential = 0;
        cellCurrent = 0;
        hydrogenProductionRate_moles = 0; # moles per second
        hydrogenProductionRate_kg = 0; # kg/s
        H2_per_hour = 0;
        predictedEfficiency = 0;
        powerUsed = 0; #MW
        powerUsedPerStack = 0;
        """
    printResults = 0;
    
    if printResults == 1:
        print("Power Used (Total): ", round(powerUsed,3), " MW")
        print("Power Used (per stack): ", round(powerUsedPerStack,3), " MW")
        print("Rated power of stack: ", round(stackRatedPower,3), " MW")
        print("Hydrogen production rate per stack: ", round(H2_per_hour/stackNumber,3)," kg/h")
        print("Hydrogen production rate total: ", round(H2_per_hour,3)," kg/h")
        print("Stack efficiency: ", round(predictedEfficiency*100,3), " %")
        print("Current per stack: ", round(cellCurrent,3), " A")
        print("Current density per stack: ", round(j/10000,3), " A/cm2")
        print("Capacity used per stack: ", round(capacityPerStack,3), " %")
        print("Voltage per stack: ", round(potential*numberOfCells,3), " V")
        print("Number of Cells per stack: ", round(numberOfCells,3), "")
    
    result = {
    "PowerUsed_MW": np.round(powerUsed, 3),
    "PowerUsedPerStack_MW": np.round(powerUsedPerStack, 3),
    "StackRatedPower_MW": np.round(stackRatedPower, 3),
    "H2_kg_h_per_stack": np.round(H2_per_hour/stackNumber, 3),
    "H2_kg_h_total": np.round(H2_per_hour, 3),
    "StackEfficiency_percent": np.round(predictedEfficiency * 100, 3),
    "Current_A_per_stack": np.round(cellCurrent, 3),
    "CurrentDensity_A_cm2": np.round(j/10000, 3),
    "CapacityUsed_percent": np.round(capacityPerStack, 3),
    "Voltage_V_per_stack": np.round(potential*numberOfCells, 3),
    "Cells_per_stack": numberOfCells
    }
    return(result)

    
    


    if mode == 'validate-polarisation':
        electrolyserTypes = ["PEM", "AEM", "Alkaline"]
        
        for electrolyserType in electrolyserTypes:
            print(f"Running simulation for: {electrolyserType}")
        
        
            if electrolyserType == "PEM":
                # Do PEM setup
                print("Using PEM parameters")
                
                # kinetic properties
                p.anode.j0 = 1;         # [A/m2]
                p.cathode.j0 = 1;       # [A/m2]
                p.anode.b = 0.03;       # [V]
                p.cathode.b = 0.04;     # [V]
                
                # electrode properties 
                p.anode.CLthickness = 10e-6; # [m]
                p.anode.specificArea = 1e+7; # [m2/m3]
                p.anode.E0 = 1.23; # [V]
                
                p.cathode.CLthickness = 10e-6; # [m]
                p.cathode.specificArea = 1e+7; # [m2/m3]
                p.cathode.E0 = 0; # [V]
                
                p.electrolyte.thickness = 50e-6; # [m]
                p.electrolyte.thickness2 = 0;    # [m] - additional layer thickness
                p.electrolyte.conductivity = 4; # [S/m]
                p.electrolyte.volFraction = 1;   # [S/m]
                
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
                # Do AEM setup
                print("Using AEM parameters")
                # kinetic properties
                p.anode.j0 = 0.0938;         # [A/m2]
                p.cathode.j0 = 1;       # [A/m2]
                p.anode.b = 0.05;       # [V]
                p.cathode.b = 0.04;     # [V]
                
                # electrode properties 
                p.anode.CLthickness = 50e-6; # [m]
                p.anode.specificArea = 2e+6; # [m2/m3]
                p.anode.E0 = 1.23; # [V]
                
                p.cathode.CLthickness = 50e-6; # [m]
                p.cathode.specificArea = 2e+6; # [m2/m3]
                p.cathode.E0 = 0; # [V]
                
                p.electrolyte.thickness = 70e-6; # [m]
                p.electrolyte.thickness2 = 0;    # [m] - additional layer thickness
                p.electrolyte.conductivity = 6; # [S/m]
                p.electrolyte.volFraction = 1;   # [S/m]
                
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
                # Do AEM setup
                print("Using AEM parameters")
                # kinetic properties
                p.anode.j0 = 0.0938;         # [A/m2]
                p.cathode.j0 = 0.695;       # [A/m2]
                p.anode.b = 0.05;       # [V]
                p.cathode.b = 0.0828;     # [V]
                
                # electrode properties 
                p.anode.CLthickness = 10e-6; # [m]
                p.anode.specificArea = 1e+6; # [m2/m3]
                p.anode.E0 = 1.23; # [V]
                
                p.cathode.CLthickness = 10e-6; # [m]
                p.cathode.specificArea = 1e+6; # [m2/m3]
                p.cathode.E0 = 0; # [V]
                
                p.electrolyte.thickness = 500e-6; # [m]
                p.electrolyte.thickness2 = 0;    # [m] - additional layer thickness
                p.electrolyte.volFraction = 0.5;   # [S/m]
                p.electrolyte.conductivity = 56*p.electrolyte.volFraction**(1.5); # [S/m]
                
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
                print("Unknown electrolyser type")
            
              
            
            # Cell Properties
            cell.area = 1;
            Asa = p.anode.CLthickness*p.anode.specificArea;
            Asc = p.cathode.CLthickness*p.cathode.specificArea;
            
            # Boundary Conditions
            phi_cathode = 0;
            phi_anode = 1.5;
            j = 1000;
            result = [];
            potentialResult = [];
            f = lambda b_a,j,j0_a,L,k,b_c,j0_c,phi_anode,E0a,Asa,Asc : b_a * np.log(j/(j0_a*Asa)) + j*L/k + b_c*np.log(j/(j0_c*Asc)) - phi_anode + E0a;
            df = lambda b_a,j,L,k,b_c : b_a/j + L/k + b_c/j;
            
            
            potential = phi_anode;
            
            
            
            
            # Set up parameters
            
            iterations = 5
            potential = 1.23  # Initial potential
            potentialResult = []
            currentDensityResult = []
            
            dx = 0.01;
            potential_values = np.arange(1.23, 3, dx)  # array of potentials
            
            for t in range(len(potential_values)):
                potential = potential_values[t]
                
                j = 0.5  # initial guess for j
                result = []
            
                for i in range(iterations):
                    phi_anode = potential  # Constant in this case
            
                    f1 = f(
                        p.anode.b, j, p.anode.j0,
                        p.electrolyte.thickness, p.electrolyte.conductivity,
                        p.cathode.b, p.cathode.j0, phi_anode, p.anode.E0, Asa, Asc
                    )
            
                    df1 = df(
                        p.anode.b, j,
                        p.electrolyte.thickness, p.electrolyte.conductivity,
                        p.cathode.b
                    )
            
                    jn = j - f1 / df1
                    j = jn
                    result.append(j)
            
                # Save the final current density and potential after convergence
                currentDensityResult.append(j)
                potentialResult.append(phi_anode)
            
            # Convert to NumPy arrays for plotting
            currentDensityResult = np.array(currentDensityResult)
            potentialResult = np.array(potentialResult)
            
        
        
            # Set marker color based on electrolyser type
            if electrolyserType == "PEM":
                markerColor = 'k'
            elif electrolyserType == "AEM":
                markerColor = 'g'
            elif electrolyserType == "Alkaline":
                markerColor = 'b'
            else:
                markerColor = 'r'  # default/fallback
            
            # Plot simulated result
            plt.plot(currentDensityResult/10000, potentialResult, color=markerColor, label = 'Model ' + electrolyserType)
            
            # Plot experimental data markers
            data_array = np.array(experimentalData)
            x_vals = data_array[:, 0]
            y_vals = data_array[:, 1]
            print("marker Color is:", markerColor)
            plt.plot(
                x_vals, y_vals,
                marker='o',
                linestyle='none',
                markerfacecolor=markerColor,
                markeredgecolor='k',  # Optional: adds contrast
                label='Exp ' + electrolyserType
                    )
            
            # Plot formatting
            plt.xlabel("Current Density (A/cm²)")
            plt.ylabel("Potential (V)")
            plt.title("Electrolyser Polarisation Curve")
            plt.grid(True)
            plt.xlim(0, 2)
            plt.ylim(1.23, 3)
            plt.legend()
            plt.tight_layout()
        plt.show()

" End of Polarisation Curve "




    






