function [x_hat] = ramsia(A, b, Zin)
% Matlab code by Huynh Van Luong, Email: huynh.luong@fau.de
% Copyright (c) 2015, Huynh Van Luong, Email: huynh.luong@fau.de
% Multimedia Communications and Signal Processing, University of Erlangen-Nuremberg.
% All rights reserved.
% 
% RAMSIA minimizes a convex problem
%   min {1/2||Ax - b||_2^2 + LAMBDA*sum(betaj*||Wj*(x - zj)||_1)}
%    x                               j 
%
% FUNCTIONS:
%   [x_hat] = ramsia(A, b, zj)
%
% INPUTS:
%   A   - m x n measurement matrix 
%   b   - m x 1 vector of observations/data 
%   Zin - n x 1 (J) vectors zj = Zin(:,j) of multiple side information 
%
% OUTPUT:
%   x_hat - n x 1 vector of recovered source 
%
%Please see the file LICENSE for the full text of the license.
% 
%     PUBLICATION: Huynh Van Luong, J. Seiler, A. Kaup, and S. Forchhammer, "Sparse Signal 
% 		Reconstruction with Multiple Side Information using Adaptive Weights for Multiview Sources," 
% 		in IEEE Int. Conf. on Image Processing 2016 (ICIP 2016), Phoenix, Arizona, USA, Sep. 2016.

%% Input optional parameters:
    t0 = tic;

    STOPPING_SPARSE_SUPPORT = 1;
    STOPPING_OBJECTIVE_VALUE = 2;
    STOPPING_SUBGRADIENT = 3;
    STOPPING_DEFAULT = STOPPING_OBJECTIVE_VALUE;
    stoppingCriterion = STOPPING_DEFAULT;
    maxIter = 5000; % maxilambdam number of iterations
    tolerance = 1e-5; % tolerance for stopping criterion.

    
%% Initializing optimization variables
    n = size(A, 2);
    t_k = 1; 
    
    invA = pinv(A);
    G = invA*A;
    L0 = max(eig(G));% default 1
    nIter = 0;
    c = invA*b;
    lambda0 = 0.5*L0*norm(c,inf);
    eta = 0.95;
    
    lambda_bar = 1e-5;
    % xk = zeros(n,1); % default
    xk = A\b; % initialied 
    lambda = lambda0;
    L = L0;

%% Input SIs
    J = size(Zin,2) + 1;
    Z = zeros(n,J);
    Z(:,1) = zeros(n,1); % Adding z0 as SI of conventional norm-1 term
    Z(:,2:end) = Zin;
%% Initialization
    timeSteps = nan(1,maxIter);
    keep_going = 1;
    nz_x = (abs(xk) > eps*10);

    Wk = zeros(n,J); % Weights on source
    Wk(:,1) = ones(size(xk,1),1); % Weights on source
    for j = 2:J
        Wk(:,j) = zeros(size(xk,1),1); % Weights on SIs
    end
    
    beta = zeros(J,1); % betas for each SI
    beta(1) = 1;
  
    f = 0.5*norm(b - A*xk)^2 + lambda_bar*sumNorm1(Wk, xk, Z); 
    xkm1 = xk;
    epsi = 0.1;%1e-1;
    epsiBeta = 1e-20;
    uk = xk;
    Wkp1 = Wk;
%% Loop step k >= 1
    % Starting loop, which is based on YALL1 package (http://yall1.blogs.rice.edu/)
    while keep_going && (nIter < maxIter)
        nIter = nIter + 1;    
        temp = G*uk - c; % gradient of f at uk  

        gk = uk - (1/L)*temp;
        xk = softMSI(gk, lambda, L, Wk, Z);

        timeSteps(nIter) = toc(t0);
        
        switch stoppingCriterion
            case STOPPING_SUBGRADIENT
                sk = L*(uk - xk) + G*(xk - uk);
                keep_going = norm(sk) > tolerance*L*max(1, norm(xk));
                
            case STOPPING_SPARSE_SUPPORT
                % compute the stopping criterion based on the change
                % of the number of non-zero components of the estimate
                nz_x_prev = nz_x;
                nz_x = (abs(xk) > eps*10);
                num_nz_x = sum(nz_x(:));
                num_changes_active = (sum(nz_x(:) ~= nz_x_prev(:)));
                if num_nz_x >= 1
                    criterionActiveSet = num_changes_active/num_nz_x;
                    keep_going = (criterionActiveSet > tolerance);
                end
                
            case STOPPING_OBJECTIVE_VALUE
                % compute the stopping criterion based on the relative
                % variation of the objective function.
                prev_f = f;
                f = 0.5*norm(b - A*xk)^2 + lambda_bar * sumNorm1(Wk, xk, Z);
                criterionObjective = abs(f - prev_f)/(prev_f);
                keep_going =  (criterionObjective > tolerance);
                
            otherwise
                error('Please define stopping criterion.');
        end   
        lambda = max(eta*lambda, lambda_bar);   

        % Weight updating
        for j = 1:J   
%             x_org = abs(xk - Z(:,j));
%             epsi = 0.1*std(x_org(x_org~=0))
            denom = (epsi + abs(xk - Z(:,j)));
            Wkp1(:,j) = 1./(denom);
            Wkp1(:,j) = Wkp1(:,j)*(n/sum(Wkp1(:,j)));
        end          
        % Beta updating
        for j = 1:J   
            denom = 0;     
            for i = 1:J
                denom = denom + (epsiBeta + sumNorm1(Wkp1(:,j), xk, Z(:,j)))/(epsiBeta + sumNorm1(Wkp1(:,i), xk, Z(:,i)));
            end
            beta(j) = 1/denom;
        end 
        % New weights
        for j = 1:J   
            Wkp1(:,j) = beta(j)*Wkp1(:,j);
        end   
        % Updating values for the next iteration
        t_kp1 = 0.5*(1 + sqrt(1 + 4*t_k*t_k));          
        ukp1 = xk + ((t_k - 1)/t_kp1)*(xk - xkm1);  
        
        % Next iteration   
        xkm1 = xk;

        Wk = Wkp1;        

        uk = ukp1;
        t_k = t_kp1;
    end
    
%% Output
    x_hat = xk;

end
%% ------------------------------------------ -----------------------------
% From this, supported subfunctions
% ------------------------------------------------------------------------

%% Soft threshoding function for multiple side information 
function y = softMSI(x, lambda, L, Wk, Z)
    A0 = -1e+20 + zeros(size(x,1),1);
    A0 = [A0, Z, A0 + 2e+20];
    [A, I] = sort(A0,2);
    W0 = [zeros(size(x,1),1), Wk, zeros(size(x,1),1)];
    W = W0;

    for i = 1:size(x,1)
        w = W0(i,:);
        W(i,:) = w(I(i,:));
    end
    y = proxMat(x, A, W, lambda, L); % Proximal function
end
%% Proximal function for softMSI
function [u] = proxMat(x, A, W, lambda, L)
    J = size(A,2) - 3;
    S = zeros(size(x,1),size(A,2) - 1);
    P = zeros(size(x,1),2*size(A,2) - 2);
    for m = 1:size(A,2) - 1
       for j = 2:J + 2
           S(:,m) = S(:,m) + W(:,j)*((-1)^(m - 2 < j - 2));       
       end
       S(:,m) = S(:,m)*(lambda/L);
       P(:,2*m - 1) = A(:,m) + S(:,m);
       P(:,2*m) = A(:,m + 1) + S(:,m);
    end
    XX = 0*P;
    
    for j = 1:2*(size(A,2) - 1)
        XX(:,j) = XX(:,j) + x;
    end
    NN = sign(P - XX);
    SNN = NN; 
    UM = SNN;
 
    for j = 2:2*(size(A,2) - 1)
        SNN(:,j) = NN(:,j-1) + NN(:,j);
    end    
    II = (SNN >= 0) .* (SNN <= 1);
    
    for j = 1:2*(size(A,2) - 1)
        if (mod(j,2) == 0)
            UM(:,j) = x - S(:,(j/2));
        else
            UM(:,j) = A(:,floor(j/2) + 1);
        end
    end      
    UM = II .* UM;
    u = sum(UM,2);
end

%% Calculating sum of norm1: sum(norm(wkj.*(xk - zj),1));
function [norm1Sum] = sumNorm1(Wk, xk, Z)
    norm1Sum = 0;
    for j = 1:size(Z,2)
        norm1Sum = norm1Sum + norm(Wk(:,j).*(xk - Z(:,j)),1);
    end
end