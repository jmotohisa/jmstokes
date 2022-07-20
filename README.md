# jmstokes

Packages for Stokes analysis and Mueller calculus. based on JM's python script

1. Stokes in ResearchReportJM
2. stokes_analysis.py in codeSnippets

# Existing related project 

[py-pol](https://py-pol.readthedocs.io/en/master/) (also at: https://bitbucket.org/optbrea/py_pol/src/master/)

## Defined functions


### Mueller matrix and calculus

+ rotT(th): Rotation matrix
    - th: angle of rotation

+ retarder(gamma): Mueller matrix of retarder with fast axis horizontal
    - gamma: retardation parameter, pi/2 for quarter waveplate, pi for half waveplate

+ rotation(m, th): Rotation of Mueller matrix
    - m: Mueller matrix
    - th: rotation angle

+ linear_polarizer(th):  Stokes parameter of the linear polarizer
    - th: angle of transmission axis

+ QWP(th): Mueller matrix of quarter wave plate
    - th: angle of rotation (0 for fast axis horizontal)

+ HWP(th): Mueller matrix of half wave plate
    - th: angle of rotation (0 for fast axis horizontal)


- MMLinPol(theta): Mueller matrix of the  Linear polarizer (with polarization angle theta)
    - theta: polarization angle in radian (theta=0 is horizontal)
- MMHWP(theta): Mueller matrix of the Half Waveplate (with polarization angle theta) (same as HWP)
- MMQWP(theta): Mueller matrix of the Quater Waveplate (with polarization angle theta) (same as QWP)

### Funcions for Stokes analysis
- CalcStokesParams(w1, w2, w3, w4, w5, w6, theta_pol=0, fnorm=False):
    Calculate stokes parameters from raw data
    - w1: lambda/2 = 0 deg
    - w2: lambda/2 = 22.5 deg
    - w3="S"+basename+"_90_none"  // lambda/2 = 45 deg
    - w4="S"+basename+"_135_none" // lambda/2 = 67.5 deg
    - w5="S"+basename+"_none_45"  // lambda/4 = 45
    - w6="S"+basename+"_none_135" // lambda/4 = 135 (-45)
    - theta_pol: angle of polarizer (1: vertical, 0: horizontal)
    - fnorm: do normalization

- ShowStokesParams2D(s0, s1, s2, s3): Display Stokes parameters for 2D matrix
    - s0, s1, s2, s3: Stokes parameters (2D matrix data)


-  SfromExEy(Ex, Ey, fnorm=False): Calculate Stokes parameters from Ex and Ey
    - Ex, Ey
	- fnorm : 
	- returns Stokes parameters


- MulMMSP(mm, s0, s1, s2, s3)
    - product of Mueller matrix and Stokes parameters (consisting of 4 waves)
	(for product of Mueller matrix, use MatrixMultiply and copy M_product to desired wave)
    - MM: muller matirx
	- returns Stokes parameters

### Analytical results of Mueller calculus

+ lin_hwp_linh(th, rhoh):
    - Input wave: Lineary polarlized light -> HWP -> horizontal linear polarizer
    - th: polarization angle of input wave
    - rhoh : rotation angle of Half Wave Plate
	
+ lin_hwp_linv(th, rhoh):
    - Input wave: Lineary polarlized light -> HWP -> vertical linear polarizer
    - th: polarization angle of input wave
    - rhoh : rotation angle of Half Wave Plate

+ circ_hwp_lin0(rhoh):
    '''
    Input wave: Circularly (R/L) polarlized light -> HWP -> vertical linear polarizer
    rhoh : rotation angle of Half Wave Plate
    '''
+ circ_hwp_lin = np.vectorize(circ_hwp_lin0, otypes=[float])

+ lin_qwp_linh(th, rhoq):
    - Input wave: Lineary polarlized light -> QWP -> horizontal linear polarizer
    - th: polarization angle of input wave
    - rhoq : rotation angle of Quarter Wave Plate

+ lin_qwp_linv(th, rhoq):
    - Input wave: Lineary polarlized light -> QWP -> vertical linear polarizer
    - th: polarization angle of input wave
    - rhoq : rotation angle of Quarter Wave Plate

+ rcirc_qwp_linh(rhoq):
    - Input wave: Right circularly polarlized light -> QWP -> horizontal linear polarizer
    - th: polarization angle of input wave
    - rhoq : rotation angle of Quarter Wave Plate

+ rcirc_qwp_linv(rhoq):
    - Input wave: Right circularly polarlized light -> QWP -> vertical linear polarizer
    - th: polarization angle of input wave
    - rhoq : rotation angle of Quarter Wave Plate

+ lcirc_qwp_linh(rhoq):
    - Input wave: Left circularly polarlized light -> QWP -> horizontal linear polarizer
    - th: polarization angle of input wave
    - rhoq : rotation angle of Quarter Wave Plate

+ lcirc_qwp_linv(rhoq):
    - Input wave: Left circularly polarlized light -> QWP -> vertical linear polarizer
    - th: polarization angle of input wave
    - rhoq : rotation angle of Quarter Wave Plate


+ lin_hwp_qwp_linh(th, rhoh, rhoq):
    - Input wave: Lineary polarlized light -> HWP -> QWP -> horizontal linear polarizer
    - th: polarization angle of input wave
    - rhoh : rotation angle of Half Wave Plate
    - rhoq : rotation angle of Quarter Wave Plate

+ linearly polarized light -> HWP -> QWP -> vertical linear polarizer


+ lin_hwp_qwp_linv(th, rhoh, rhoq):
    - Input wave: Lineary polarlized light -> HWP -> QWP -> vertical linear polarizer
    - th: polarization angle of input wave
    - rhoh : rotation angle of Half Wave Plate
    - rhoq : rotation angle of Quarter Wave Plate

+ rcirc_hwp_qwp_linh(rhoh, rhoq):
    - Input wave: R (L) circularly polarlized light -> HWP -> QWP -> horizontal (vertical) linear polarizer
    - rhoh : rotation angle of Half Wave Plate
    - rhoq : rotation angle of Quarter Wave Plate


+ lcirc_hwp_qwp_linh(rhoh, rhoq):

    - Input wave: L (R) circularly polarlized light -> HWP -> QWP -> horizontal (vertical) linear polarizer
    - rhoh : rotation angle of Half Wave Plate
    - rhoq : rotation angle of Quarter Wave Plate


- linearly polarized light -> HWP -> QWP -> vertical linear polarizer
+ lin_qwp_hwp_linh(th, rhoq, rhoh):
    - Input wave: Lineary polarlized light -> QWP -> HWP -> horizontal linear polarizer
    - th: polarization angle of input wave
    - rhoh : rotation angle of Half Wave Plate
    - rhoq : rotation angle of Quarter Wave Plate

+ rcirc_qwp_hwp_linh(rhoq, rhoh):
    - Input wave: R (L) circularly polarlized light -> QWP -> HWP -> horizontal (vertical) linear polarizer
    - rhoh : rotation angle of Half Wave Plate
    - rhoq : rotation angle of Quarter Wave Plate

