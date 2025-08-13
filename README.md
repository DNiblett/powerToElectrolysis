# powerToElectrolysis
Code to take power from renewable energy (e.g. wind) and convert it to hydrogen production rate (kg/hr) using a generalised electrolyser model for PEM, Alkaline and AEM. Time supplied in seconds and power in MW. Example file taken from: Niblett, Daniel; Yeter, Baran; Mamlouk, Mohamed (2023). Wind Speed & Power Generated Dataset For Floating Offshore 7 MW and 15 MW Turbine. Newcastle University. Dataset. https://doi.org/10.25405/data.ncl.24516718.v1

1. runPowerToElectrolysis.py - this reads the time series data from the .csv file and then runs "electrolyser_full.py"
2. electrolyser_full.py - using the electrolyser type (PEM,Alkaline,AEM), number of stacks, stack rated capacity (MW), power input (MW) and minimum load (%) this solves for the current at each time. Specific electrolyser materials, catalysts, thicknesses can be changed within the script.
3. Stack rated capacity is used along with cell area (fixed inside each electrolyte type) to calculate number of cells.
4. Minimum load is used as post-processing to set all parameters to 0 when the power applied to each stack is less than the minimum load.
5. electrolyserPolarisation.py - this code produces the polarisation curves for the three electrolysers, with parameters roughly chosen to fit some experimental data in literature.
