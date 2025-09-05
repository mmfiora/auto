# -*- coding: utf-8 -*-
"""
Created on Fri Sep  5 16:45:01 2025

@author: Maria
"""

import pandas as pd
import matplotlib.pyplot as plt

df=pd.read_csv("activity_normalized.csv")

species_counts=df["species"].value_counts()

plt.figure(figsize=(14, 8))
species_counts.plot(kind="bar")


plt.xlabel("Species", fontsize=12)
plt.ylabel("Frequency", fontsize=12)
plt.xticks(rotation=90)  # etiquetas verticales
plt.tight_layout()
plt.show()