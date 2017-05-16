# RAMSIA
Reconstruction Algorithm with Multiple Side Information using Adaptive Weights (RAMSIA)

	Version 1.1, January 15, 2016
	Implementations by Huynh Van Luong, Email: huynh.luong@fau.de,
	Multimedia Communications and Signal Processing, University of Erlangen-Nuremberg.

`Please see` [LICENSE](https://github.com/huynhlvd/ramsia/blob/master/LICENSE.md) `for the full text of the license.`

    PUBLICATION: Huynh Van Luong, J. Seiler, A. Kaup, and S. Forchhammer, "Sparse Signal 
		Reconstruction with Multiple Side Information using Adaptive Weights for Multiview Sources," 
		in IEEE Int. Conf. on Image Processing 2016 (ICIP 2016), Phoenix, Arizona, USA, Sep. 2016.

  **_Solving the _n-l1_ minimization problem:_**
  
<img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\min_{\boldsymbol{x}}\Big\{H(\boldsymbol{x})=&space;\frac{1}{2}\|\mathbf{\Phi}\boldsymbol{x}-\boldsymbol{y}\|^{2}_{2}&space;&plus;&space;\lambda&space;\sum\limits_{j=0}^{J}\beta_j\|\mathbf{W}_{j}(\boldsymbol{x}-\boldsymbol{z}_{j})\|_{1}\Big\}" title="\min_{\boldsymbol{x}}\Big\{H(\boldsymbol{x})= \frac{1}{2}\|\mathbf{\Phi}\boldsymbol{x}-\boldsymbol{y}\|^{2}_{2} + \lambda \sum\limits_{j=0}^{J}\beta_j\|\mathbf{W}_{j}(\boldsymbol{x}-\boldsymbol{z}_{j})\|_{1}\Big\}" />

Inputs:
- <img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\boldsymbol{y}\in&space;\mathbb{R}^{m}" title="\boldsymbol{y}_{t}\in \mathbb{R}^{m}" />: A vector of observations/data <br /> 
- <img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\mathbf{\Phi}\in&space;\mathbb{R}^{m\times&space;n}" title="\mathbf{\Phi}\in \mathbb{R}^{m\times n}" />: A measurement matrix <br />
- <img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\boldsymbol{z}_{j}\in&space;\mathbb{R}^{n}" title="\boldsymbol{y}_{t}\in \mathbb{R}^{m}" />: Multiple side information signals <br />

Outputs:
- <img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\boldsymbol{\hat{x}}\in\mathbb{R}^{n}" title="\boldsymbol{x}_{t},\boldsymbol{v}_{t}\in\mathbb{R}^{n}" />: The recovered source

**_Source code files:_**  
 - [ramsia.m](https://github.com/huynhlvd/ramsia/blob/master/ramsia.m): The function for RAMSI
 - [usageDemo_ramsia.m](https://github.com/huynhlvd/ramsia/blob/master/usageDemo_ramsia.m): One demo to run RAMSI
