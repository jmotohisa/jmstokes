#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import muellerCalculus as mc

polh0 = np.matrix([[1, 1, 0, 0], [1, 1, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]]) / \
    2.  # Mueller matrix of polarizer with transmission axis horizontal
# Mueller matrix of half wave plate with fast axis horizontal
hwp0 = mc.retarder(math.pi)
# Mueller matrix of quarter wave plate with fast axis horizontal
qwp0 = mc.retarder(math.pi/2)
# Stokes parameter for horizontally linearly polarized light
s_linpol0 = [1, 1, 0, 0]
s_rpol0 = [1, 0, 0, 1]  # Stokes parameter for right ciruclarly polarized light

sin = [1, -0.1, 0, 0.08]
th_list = np.linspace(0, math.pi, endpoint=True)

res1_list = np.dot(np.dot(mc.linear_polarizer(math.pi/2),
                          mc.rotation(hwp0, th_list)), sin)[0, 0]
res2_list = np.dot(np.dot(mc.linear_polarizer(math.pi/2),
                          mc.rotation(qwp0, th_list)), sin)[0, 0]


fig, ax = plt.subplots()
# ax.xaxis.set_minor_locator(MultipleLocator(5))
ax.plot(th_list*180/math.pi, res1_list, label='HWP', color='magenta')
ax.plot(th_list*180/math.pi, res2_list, label='QWP', color='green')
ax.legend(loc='best', fontsize=18)
ax.set_title('HWP/QWP+linear horizontal polarizer', fontsize=18)
ax.set_xlabel('HWP/QWP angle (degree)', fontsize=18)
ax.set_ylabel('intensity', fontsize=18)
ax.tick_params(labelsize=18)
ax.set_xticks([0, 45, 90, 135, 180])
# ax.xaxis.set_minor_locator(MultipleLocator(15))
ax.set_xlim(0, 180)
ax.set_ylim(0, 1)
plt.show()
