Reconstruction Algorithm with Multiple Side Information using Adaptive Weights(RAMSIA)
Version 1.1, January 15, 2016
Implementations by Huynh Van Luong, Email: huynh.luong@fau.de

Please see the file LICENSE for the full text of the license.
-----------------------------------------------------------------------------------------
The main files are:
    */ramsia.m: The function for RAMSIA
    */usageDemo.m: One demo to run RAMSIA
Some simple changes may be made to each of these files to get them to run on your machine. 

PUBLICATION: Huynh Van Luong, J. Seiler, A. Kaup, and S. Forchhammer, "Sparse Signal 
Reconstruction with Multiple Side Information using Adaptive Weights for Multiview Sources," 
in IEEE Int. Conf. on Image Processing 2016 (ICIP 2016), Phoenix, Arizona, USA, Sep. 2016.

Solving the n-l1 minimization problem:
	minimize    \sum betaj*||Wj*(x - zj)||_1 subject to  Ax = b 
       x	       j
    where b: m x 1,  A: m x n,  x, zj: n x 1, Wj = diag(n), betaj > 0, wij (\in Wj) > 0.
        Implementations by Huynh Van Luong, Email: huynh.luong@fau.de
        Multimedia Communications and Signal Processing, University of Erlangen-Nuremberg.
        
        Please see the file LICENSE for the full text of the license.
-----------------------------------------------------------------------

**I. Reconstruction Algorithm with Multiple Side Information (RAMSI)**

    PUBLICATION: H. V. Luong, J. Seiler, A. Kaup, and S. Forchhammer, "A Reconstruction Algorithm with 
    	Multiple Side Information for Distributed Compression of Sparse Sources," in Data Compression 
    	Conference 2016 (DCC 2016), Snowbird, Utah, Apr. 2016.

  **_Solving the _n-l1_ minimization problem:_**
  
<img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\min_{\boldsymbol{x}}\Big\{H(\boldsymbol{x})=&space;\frac{1}{2}\|\mathbf{\Phi}\boldsymbol{x}-\boldsymbol{y}\|^{2}_{2}&space;&plus;&space;\lambda&space;\sum\limits_{j=0}^{J}\|\mathbf{W}_{j}(\boldsymbol{x}-\boldsymbol{z}_{j})\|_{1}\Big\}" title="\min_{\boldsymbol{x}}\Big\{H(\boldsymbol{x})= \frac{1}{2}\|\mathbf{\Phi}\boldsymbol{x}-\boldsymbol{y}\|^{2}_{2} + \lambda \sum\limits_{j=0}^{J}\|\mathbf{W}_{j}(\boldsymbol{x}-\boldsymbol{z}_{j})\|_{1}\Big\}"  /> (1)

Inputs:
- <img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\boldsymbol{y}\in&space;\mathbb{R}^{m}" title="\boldsymbol{y}_{t}\in \mathbb{R}^{m}" />: A vector of observations/data <br /> 
- <img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\mathbf{\Phi}\in&space;\mathbb{R}^{m\times&space;n}" title="\mathbf{\Phi}\in \mathbb{R}^{m\times n}" />: A measurement matrix <br />
- <img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\boldsymbol{z}_{j}\in&space;\mathbb{R}^{n}" title="\boldsymbol{y}_{t}\in \mathbb{R}^{m}" />: Multiple side information signals <br />

Outputs:
- <img src="https://latex.codecogs.com/svg.latex?\dpi{150}&space;\boldsymbol{x}\in\mathbb{R}^{n}" title="\boldsymbol{x}_{t},\boldsymbol{v}_{t}\in\mathbb{R}^{n}" />: Estimates of foreground and background

**_Source code files:_**  
 - `ramsi.m`: The function for RAMSI
 - `usageDemo_ramsi.m`: One demo to run RAMSI
