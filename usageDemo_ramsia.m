% Demo of using ramsia.m for RAMSIA
% Matlab code by Huynh Van Luong, Email: huynh.luong@fau.de
% Jan. 15, 2016
% 
%% Initialization 
% A  - m x n measurement matrix 
n = 1000;
m = 130;
A = randn(m,n);
% Supports of ||x - zj||_0 = sj
s0 = 128;
sj = 64;
s1 = sj; 
s2 = sj;
s3 = sj;
S = [s1 s2 s3]; 
J = size(S,2);
%% Generating x
ranDat = [randn(s0,1); zeros(n-s0,1)];
perm = randperm(n);
x = ranDat(perm);
%% Generating zj = Z(:,j)
Z = zeros(n,J);
for j = 1:J
    ran = [randn(S(j),1); zeros(n - S(j),1)];
    ran = ran(perm);
    Z(:,j) = x + ran; % Positions of non-zeros of x and zj are coincided 
end

%% Running RAMSIA
% Input observation b with m measurements
b = A*x;
x_hat = ramsia(A, b, Z);
er(i) = norm(x_hat - x,2)/(norm(x,2));
fprintf('Recovered error: %4.8f \n', norm(x_hat - x,2)/(norm(x,2)));

