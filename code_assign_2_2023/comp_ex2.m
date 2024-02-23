%% 4 RQ Factorization and Computation of K

%%% Computer Exercise 2 %%%
clc;clearvars;close all;

load for_comp_e2.mat
load compEx1data.mat

% Task:
% Compute K for one camera in each of the solutions you obtained in computer exercise 1.
% Choosing camera P(1)
% Submit the K matrices. (Make sure that element K(3,3) = 1 by division; K = K./K(3, 3).)


% Re-generate camera matrices from Ex1
P1 = P{1};
PT1 = P{1}*inv(T1);
PT2 = P{1}*inv(T2);   

% Calculate K
% RQ factorization
% r upper tri, q unit matrix
[K1,q1]     = rq(P1);
[K1T1,qT1]   = rq(PT1);
[K1T2,qT2]   = rq(PT2);

% Make sure K(3,3) = 1
K1  = K1./K1(3,3);
K1T1 = K1T1./K1T1(3,3);
K1T2 = K1T2./K1T2(3,3);


% Do they represent the same transformation?
test1 = K1 == K1T1  % not the same
test2 = K1T1 == K1T2 % not the same
test3 = K1 == K1T2  % the same!


save('K_matrices.mat', 'K1', 'K1T1', 'K1T2');


%%% Answer %%%
% In the case of camera matrix P(1) and P(1) with T2 transformation, yes
% they do represent the same transformation. This is not true for all
% combinations of P(1) with T1 transformation. This is not due to a scale
% factor not being correct. 


