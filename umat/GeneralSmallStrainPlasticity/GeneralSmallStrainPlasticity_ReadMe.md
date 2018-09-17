# Rate Independent Plasticity
## Base models
* chaboche
* delobelle
* ohnowang
## Rate dependent models (which can be combined with the above):
* norate: Rate independent version of base model
* norton: Norton overstress
* cowsym: Cowper-Symonds overstress
* delobelle: Delobelle overstress


## Model description
The number of state variables are: 2 + nvaradd + 6*nback, where nvaradd=-1 for rate dependent models
The material parameters are different for each base model (additional parameters for rate dependent given in curly braces {}):
chaboche:   E, nu, Y0, Hiso, invYiso, {tstar, nexp,} Hk1, invYk1, Hk2, ...             (5+2*nback)
delobelle:  E, nu, Y0, Hiso, invYiso, delta, {tstar, nexp,} Hk1, invYk1, Hk2, ...      (6+2*nback)
ohnowang:   E, nu, Y0, Hiso, invYiso, delta, {tstar, nexp,} Hk1, invYk1, mk1, Hk2, ... (6+3*nback)

Description of the material parameters.
E:          Modulus of elasticity
nu:         Poissons ratio
Y0:         Initial yield limit
Hiso:       Isotropic hardening modulus
invYiso:    The inverse of the isotropic hardening saturation stress
delta:      Weighting between Armstrong-Frederick and Burlet-Cailletaud kinematic saturation for all back stresses
tstar:      Relaxation time
nexp:       Exponent in overstress function
Hki:        Kinematic hardening modulus nr. i
invYki:     The inverse of the kinematic hardening saturation stress nr. i
mki:        The exponent in the Ohno-Wang model for backstress nr. i

## Building the linked libraries
To build the models and create a linked library do the following
### Windows with Visual Studio
In command window:
* <mkdir build>
* <cd build>
* Windows with Visual Studio: <cmake -G "Visual Studio 11 2012 Win64" ../src> (Replace with appropriate version of VS)
* Linux: <cmake ../src>
* <cmake -D basemod=<Base Model> -D ratemod=<Rate Dependent Model> .
* <cmake --build . --config RELEASE
* Move the .dll (Windows) or .so (Linux) file for the model that you've created (found in the build/Release (Windows) or build (Linux) folder) to your material model folder (e.g. the same as the matmodfit program)
* Windows: If you need debugging capabilities, also copy the .pdb file with the same name (in which case the config DEBUG should have been used instead)
* If you want to build all the models, the script build_all can be used:
* Windows: <python bld_all_models.py "-G \"Visual Studio 11 2012 Win64\"">
* Linux: <python bld_all_models.py>

## To make subroutines for Abaqus
* Rename the GeneralSmallStrain.f90 to GeneralSmallStrain.for (.f on linux)
* Uncomment the necessary includes in the GeneralSmallStrain.f90 file
* Comment away the !DEC EXPORT line (add an additional exclamation mark: !)
* run: abaqus make library=GeneralSmallStrain.for (or .f)
* Move the generated files and reset the changes made above (should not be part of any commit)