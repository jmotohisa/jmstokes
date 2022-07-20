#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from muellerCalculus import *

polh0 = np.matrix([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]) / \
    2.  # Mueller matrix of polarizer with transmission axis horizontal
# Mueller matrix of half wave plate with fast axis horizontal
hwp0 = retarder(math.pi)
# Mueller matrix of quarter wave plate with fast axis horizontal
qwp0 = retarder(math.pi/2)
# Stokes parameter for horizontally linearly polarized light
s_linpol0 = [1, 1, 0, 0]
# Stokes parameter for right ciruclarly polarized light
s_rpol0 = [1, 0, 0, 1]

phi_m = np.linspace(0, np.pi, 100)
phi_p = np.linspace(0, np.pi, 100)
X, Y = np.meshgrid(phi_p, phi_m)

Z = np.dot(np.dot(np.dot(linear_polarizer(0), rotation(qwp0, X)),
                  rotation(hwp0, Y)), s_linpol0)[0, 0]
fig, ax = plt.subplots(figsize=(5, 4.5))
p = ax.pcolor(X/(np.pi)*180, Y/(np.pi)*180, Z,
              cmap=matplotlib.cm.RdBu, vmin=0, vmax=1)
cb = fig.colorbar(p, ax=ax)
ax.set_title(
    'Input: horizontal linear polarization\nHWP+QWP+linear horizontal pol.', fontsize=18)
ax.set_xlabel('QWP angle (degree)', fontsize=18)
ax.set_ylabel('HWP angle (degree)', fontsize=18)
ax.tick_params(labelsize=18)
ax.set_xticks([0, 45, 90, 135, 180])
ax.set_xlim(0, 180)
ax.set_yticks([0, 45, 90, 135, 180])
ax.set_ylim(0, 180)
plt.show()

Z = np.dot(np.dot(np.dot(linear_polarizer(0), rotation(hwp0, Y)),
                  rotation(qwp0, X)), s_rpol0)[0, 0]
fig, ax = plt.subplots(figsize=(5, 4.5))
p = ax.pcolor(X/(np.pi)*180, Y/(np.pi)*180, Z,
              cmap=matplotlib.cm.RdBu, vmin=0, vmax=1)
cb = fig.colorbar(p, ax=ax)
ax.set_title(
    'Input: horizontal linear polarization\nHWP+QWP+linear horizontal pol.', fontsize=18)
ax.set_xlabel('QWP angle (degree)', fontsize=18)
ax.set_ylabel('HWP angle (degree)', fontsize=18)
ax.tick_params(labelsize=18)
ax.set_xticks([0, 45, 90, 135, 180])
ax.set_xlim(0, 180)
ax.set_yticks([0, 45, 90, 135, 180])
ax.set_ylim(0, 180)
plt.show()
