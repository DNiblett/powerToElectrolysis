#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 14:22:52 2025

@author: Dr Daniel Niblett, Senior Research Associate, Newcastle University

Code which takes .csv input of power in MW from turbine, solves for electrolyser
current based on power input and design of electrolyser.

"""
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from electrolyser_full import runElectrolyserFull

# ---------------------------
# Load CSV (with column names)
# ---------------------------
data = pd.read_csv("Data_7MW_Turbine.csv", header=None, names=["time_s", "wind_m_s", "power_MW"])

# Prefer using columns directly (clearer, avoids index mixups)
time_s   = data["time_s"].to_numpy()
wind_m_s = data["wind_m_s"].to_numpy()
powerInput = data["power_MW"].to_numpy()

# ---------------------------
# Model settings
# ---------------------------
stackNumber   = 5           # number of stacks
stackCapacity = 2           # MW per stack
minLoad       = 10          # minimum load in %
electrolyserType = "Alkaline"    # "PEM", "AEM", "Alkaline" (case-sensitive)
runMode = "transient"

# ---------------------------
# Run electrolyser model once over the whole input vector
# (assuming runElectrolyserFull accepts vector powerInput)
# ---------------------------
out = runElectrolyserFull(powerInput, stackNumber, stackCapacity, minLoad, electrolyserType, runMode)

# ---------------------------
# Clean NaNs and align scalars to array length
# ---------------------------
N = len(powerInput)

def expand_scalar_to_array(val, n):
    if np.isscalar(val):
        return np.full(n, val)
    val = np.asarray(val)
    if val.shape[0] != n:
        raise ValueError(f"Output array length {val.shape[0]} doesn't match input length {n}.")
    return val

def nan_to_zero(val):
    # Works for scalars and arrays
    if val is None:
        return 0
    arr = np.asarray(val)
    if arr.dtype.kind in "fc":  # float/complex
        arr = np.nan_to_num(arr, nan=0.0, posinf=0.0, neginf=0.0)
    return arr

out_aligned = {}
for k, v in out.items():
    v_expanded = expand_scalar_to_array(v, N)
    out_aligned[k] = nan_to_zero(v_expanded)

# ---------------------------
# Build a single tidy DataFrame
# ---------------------------
df = pd.DataFrame({
    "time_s": time_s,
    "wind_m_s": wind_m_s,
    **out_aligned
})

# ---------------------------
# Mask rows where CapacityUsed_percent < minLoad
# Keep time_s and wind_m_s unchanged; zero out computed outputs
# ---------------------------
cols_to_zero = [c for c in df.columns if c not in ("time_s", "wind_m_s")]
df.loc[df["CapacityUsed_percent"] < minLoad, cols_to_zero] = 0

# ---------------------------
# Basic calculations and plots
# ---------------------------
h2_rate = df["H2_kg_h_total"].to_numpy()               # kg/h
dt = np.diff(time_s, prepend=time_s[0])                # seconds
accumulated_H2 = (h2_rate / 3600.0) * dt               # kg produced during each dt
total_H2 = np.cumsum(accumulated_H2)                   # kg
final_H2 = float(total_H2[-1])
print("Total Hydrogen Produced (kg):", final_H2)

hours = (time_s - time_s.min()) / 3600.0

plt.figure()
plt.plot(hours, powerInput, linewidth=0.5)
plt.xlabel("Time (hours)")
plt.ylabel("Power Input (MW)")
plt.title("Power Input")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(hours, total_H2, linewidth=0.5)
plt.xlabel("Time (hours)")
plt.ylabel("Total H2 produced (kg)")
plt.title("H2 Production Total")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(hours, h2_rate, linewidth=0.5)
plt.xlabel("Time (hours)")
plt.ylabel("H2 production rate (kg/h)")
plt.title("H2 Production Rate")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(hours, df["StackEfficiency_percent"], linewidth=0.5)
plt.xlabel("Time (hours)")
plt.ylabel("Efficiency (%)")
plt.title("Stack Efficiency")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(hours, df["Voltage_V_per_stack"], linewidth=0.5)
plt.xlabel("Time (hours)")
plt.ylabel("Stack Voltage (V)")
plt.title("Stack Voltage")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(hours, df["CapacityUsed_percent"], linewidth=0.5)
plt.xlabel("Time (hours)")
plt.ylabel("Capacity Per Stack (%)")
plt.title("Capacity Used")
plt.grid(True)
plt.tight_layout()
plt.show()

plt.figure()
plt.plot(hours, df["Current_A_per_stack"], linewidth=0.5)
plt.xlabel("Time (hours)")
plt.ylabel("Current per stack (A)")
plt.title("Current per Stack")
plt.grid(True)
plt.tight_layout()
plt.show()

# ---------------------------
# Save to CSV
# ---------------------------
fileName = "Electrolyser_Output_" + electrolyserType + ".csv"
df.to_csv(fileName, index=False)
print("Saved",fileName)
