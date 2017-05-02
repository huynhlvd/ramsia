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
