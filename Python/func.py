# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
import matplotlib.pyplot as plt

a = 2

x = np.linspace(-7,7,1000)
func = np.exp(-a*np.abs(x))
plt.plot(x,func)
plt.axis([-6,6, 0,1])
plt.title("Plot of function to integrate in the x dimension")
plt.xlabel("x")
plt.ylabel("$e^{-\\alpha \\sqrt{x^2_i}}$", fontsize = 14)
plt.savefig("/Users/andreas/Computational Physics/Fys4150/Project 3/Python/figure1")