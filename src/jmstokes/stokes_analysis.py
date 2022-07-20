#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt


def CalcStokesParams(w1, w2, w3, w4, w5, w6, theta_pol=0, fnorm=False):
    """
    Calculate stokes parameters from raw data
        w1: lambda/2 = 0 deg
        w2: lambda/2 = 22.5 deg
        w3="S"+basename+"_90_none"  // lambda/2 = 45 deg
        w4="S"+basename+"_135_none" // lambda/2 = 67.5 deg
        w5="S"+basename+"_none_45"  // lambda/4 = 45
        w6="S"+basename+"_none_135" // lambda/4 = 135 (-45)
    theta_pol: angle of polarizer (1: vertical, 0: horizontal)
    fnorm: do normalization
    """
    s0 = np.empty_like(w1)
    s1 = np.empty_like(w1)
    s2 = np.empty_like(w1)
    s3 = np.empty_like(w1)

    if(theta_pol == 1):
        s0 = (w3+w1)
        s1 = (w3-w1)
        s2 = (w4-w2)
        s3 = (w5-w6)
    else:
        s0 = (w1+w3)
        s1 = (w1-w3)
        s2 = (w2-w4)
        s3 = (w6-w5)

    if(fnorm):
        s1 = s1/s0
        s2 = s2/s0
        s3 = s3/s0

    return s0, s1, s2, s3


def ShowStokesParams2D(s0, s1, s2, s3):
    """
    Display Stokes parameters for 2D matrix
    """

    fig, axs = plt.subplots(4, 1)

    ax1 = axs[0]
    pcm1 = ax1.pcolor(X, Y, s0, cmap='hsv', shading='auto')
    fig.colorbar(pcm1, ax=ax1, orientation="vertical")

    ax2 = axs[1]
    pcm2 = ax2.pcolor(X, Y, s1, cmap='hsv', shading='auto')
    fig.colorbar(pcm2, ax=ax2, orientation="vertical")

    ax3 = axs[2]
    pcm3 = ax3.pcolor(X, Y, s2, cmap='hsv', shading='auto')
    fig.colorbar(pcm3, ax=ax3, orientation="vertical")

    ax4 = axs[3]
    pcm4 = ax4.pcolor(X, Y, s3, cmap='hsv', shading='auto')
    fig.colorbar(pcm4, ax=ax4, orientation="vertical")

    plt.show()

    return


def SfromExEy(Ex, Ey, fnorm=False):
    """
    Calculate Stokes parameters from Ex and Ey
    """
    s0 = np.empty_like(Ex)
    s1 = np.empty_like(Ex)
    s2 = np.empty_like(Ex)
    s3 = np.empty_like(Ex)

    s0 = Ex**2+Ey**2
    s1 = Ex**2-Ey**2
    s2 = 2*Ex*Ey
    s3 = 0

    if(fnorm):
        s1 /= s0
        s2 /= s0

    return s0, s1, s2, s3


def MMLinPol(theta):
    """
    Mueller matrix
    Linear polarizer (with polarization angle theta)
    theta: polarization angle in radian (theta=0 is horizontal)
    """
    mm = np.zeros((4, 4))
    th2 = theta*2

    mm[0, 0] = 1
    mm[1, 0] = np.cos(th2)
    mm[2, 0] = np.sin(th2)
    mm[3, 0] = 0
    mm[0, 1] = np.cos(th2)
    mm[1, 1] = np.cos(th2)**2
    mm[2, 1] = np.cos(th2)*np.sin(th2)
    mm[3, 1] = 0
    mm[0, 2] = np.sin(th2)
    mm[1, 2] = np.cos(th2)*np.sin(th2)
    mm[2, 2] = np.sin(th2)**2
    mm[3, 2] = 0
    mm[0, 3] = 0
    mm[1, 3] = 0
    mm[2, 3] = 0
    mm[3, 3] = 0

    mm = mm/2

    return mm


def MMQWP(theta):
    """
    Mueller matrix
    Quater Waveplate (with polarization angle theta)
    """
    mm = np.zeros((4, 4))
    th2 = theta*2

    mm[0, 0] = 1
    mm[1, 0] = 0
    mm[2, 0] = 0
    mm[3, 0] = 0
    mm[0, 1] = 0
    mm[1, 1] = np.cos(th2)**2
    mm[2, 1] = np.cos(th2)*np.sin(th2)
    mm[3, 1] = np.sin(th2)
    mm[0, 2] = 0
    mm[1, 2] = np.cos(th2)*np.sin(th2)
    mm[2, 2] = np.sin(th2)**2
    mm[3, 2] = -np.cos(th2)
    mm[0, 3] = 0
    mm[1, 3] = -np.sin(th2)
    mm[2, 3] = np.cos(th2)
    mm[3, 3] = 0

    return mm


def MMHWP(theta):
    """
    Mueller matrix
    Half Waveplate (with polarization angle theta)
    """
    th4 = theta*4
    mm = np.zeros((4, 4))

    mm[0, 0] = 1
    mm[1, 0] = 0
    mm[2, 0] = 0
    mm[3, 0] = 0
    mm[0, 1] = 0
    mm[1, 1] = np.cos(th4)
    mm[2, 1] = np.sin(th4)
    mm[3, 1] = 0
    mm[0, 2] = 0
    mm[1, 2] = np.sin(th4)
    mm[2, 2] = -np.cos(th4)
    mm[3, 2] = 0
    mm[0, 3] = 0
    mm[1, 3] = 0
    mm[2, 3] = 0
    mm[3, 3] = -1

    return mm


def MulMMSP(mm, s0, s1, s2, s3):
    """
    product of Mueller matrix and Stokes parameters (consisting of 4 waves)
    (for product of Mueller matrix, use MatrixMultiply and copy M_product to desired wave)
    MM: muller matirx
    """
    s0_dest = m[0, 0]*s0 + m[0, 1]*s1 + m[0, 2]*s2 + m[0, 3]*s3
    s1_dest = m[1, 0]*s0 + m[1, 1]*s1 + m[1, 2]*s2 + m[1, 3]*s3
    s2_dest = m[2, 0]*s0 + m[2, 1]*s1 + m[2, 2]*s2 + m[2, 3]*s3
    s3_dest = m[3, 0]*s0 + m[3, 1]*s1 + m[3, 2]*s2 + m[3, 3]*s3

    return s0_dest, s1_dest, s2_dest, s3_dest
