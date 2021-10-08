# Mancinelli_2021_PNAS
This repo contains the matlab files referenced in the publication from
Mancinelli et al. published 2021 in PNAS.
(*Link will be included as soon as the paper is published*)

## Installation
clone this repo or download files in order to use it

## Description
This repo contains a reductionist Monte Carlo model based on the
Gillespie algorithm to simulate curvature-dependent protein-membrane
interactions as well as actin filament outgrowth in neuronal dendrites.

**The model was implemented for various geometries in 2D and 3D:**
* Mancinelli_2021_PNAS_2DLine.m simulates the longitudinal section of a
neuronal dendrite
* Mancinelli_2021_PNAS_2DCircle.m simulates the cross section of a
neuronal dendrite
* Mancinelli_2021_PNAS_3DCylinder.m simulates a cylindrical section of a
neuronal dendrite

For more details please see the supplementary materials.
(*Link will be included as soon as the paper is published*)

## Usage
To execute the simulation with default parameters, just call the script
```matlab
Mancinelli_2021_PNAS_2DLine
```
Alternativly you can hit the run button. In order to use custom input
parameters, call the script and add the parameter and the desired
parameter value in the brackets.
```matlab
Mancinelli_2021_PNAS_2DLine('membrane rigidity',20)
```
Alternativly you can change the default settings directly in the source
code.
Some parameters apply for all three, some are specific to one model.
Description, units and the physiological range of the parameters are
available in table 2-5 in the supplementary materials. (*Link will be
added as soon as the paper is published*)
The parameter settings used in Maninelli et al. can be found in
*Mancinelli_2021_PNAS_ParameterSettings.pdf*.


#### general parameters:
##### membrane parameters

| parameter name  | default value | additional description |
| --------------- |:-------------:|-------------|
| 'D'             | 3             | membrane diffusion constant |
|'dh'|1|diffusion step size|
| 'membrane stretch modulus'| 0     |membrane stretch modulus|
| 'membrane rigidity' | 10     |membrane rigidity|
| 'membrane c0' | 0 | membrane intrinsic curvature |

##### I-BAR parameters

| parameter name  | default value | additional description |
| --------------- |:-------------:|-------------|
|'Ibar k1'|0.11|adsorption constant|
|'Ibar k2'|36|desorption constant|
|'Ibar c0' | -0.05 |intrinsic curvture |
|'Ibar rigidity'|10|rigidity |
|'Ibar area'|35|protein cross section|
|'Ibar standard state concentration'|1|standard state concentration|
|'Ibar mu'|1|membrane affinity|
|'Ibar mu0'|0|standard state chemical potential|
|'Ibar start saturation'|0|start saturation|

##### N-BAR parameters

| parameter name  | default value | additional description |
| --------------- |:-------------:|-------------|
|'Nbar k1'|0.11|adsorption constant|
|'Nbar k2'|36|desorption constant|
|'Nbar c0' | 0.05 |intrinsic curvture |
|'Nbar rigidity'|10|rigidity |
|'Nbar area'|35|protein cross section|
|'Nbar standard state concentration'|1|standard state concentration|
|'Nbar mu'|1|membrane affinity|
|'Nbar mu0'|0|standard state chemical potential|
|'Nbar start saturation'|0|start saturation|

##### actin parameters

| parameter name  | default value | additional description |
| --------------- |:-------------:|-------------|
|'actin k_on'|11.6|polymerization rate constant|
|'actin k_off'|1.4|depolymerization rate constant|
|'actin concentration'|10|cytosolic G-actin concentration|
|'actin start position'|-50|filaments' initial distance to membrane|

##### numeric parameters
| parameter name  | default value | additional description |
| --------------- |:-------------:|-------------|
|'tend'|60|total simulation time|
|'plot on'|true|plot current state|
|'plot everyt'|2|plot current state every x seconds|
|'save plot'|false|save current state|
|'save everyt'|2|save current state every x seconds|

#### Mancinelli_2021_PNAS_2DLine specific parameters:
| parameter name  | default value | additional description |
| --------------- |:-------------:|-------------|
|'starting shape'|'flat'|initial shape of the membrane|
|'L'|1000|length of simulated domain|
|'dL'|15|membrane patch length|
|'membrane width'|50|membrane patch width|

#### Mancinelli_2021_PNAS_2DCircle specific parameters:
| parameter name  | default value | additional description |
| --------------- |:-------------:|-------------|
|'R'|1000|circle radius|
|'dR'|15|membrane patch arc length|
|'membrane width'|50|membrane patch width|

#### Mancinelli_2021_PNAS_3DCylinder specific parameters:
no input parser.

## Energy Landscape
The *EnergyLandscape.mlapp* is a GUI that shall help to predict local energy minima for a specific parameter combination. To start the GUI call the file from console or hit the run button.

## Abstract Mancinelli et al. 2021
Abstract will be included as soon as the paper is published

## License
[MIT](https://choosealicense.com/licenses/mit/)
