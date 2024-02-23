%% 6. Feature Extraction and Matching using SIFT
clear;clc;close all;

%%% Computer Exercise 4 %%%

im1 = imread('cube1.JPG');
figure(1)
imagesc(im1)
im2 = imread('cube2.JPG');
figure(2)
imagesc(im2)


[f1 d1] = vl_sift( single(rgb2gray(im1)), 'PeakThresh', 1);
[f2 d2] = vl_sift( single(rgb2gray(im2)), 'PeakThresh', 1);


figure(1)
vl_plotframe(f1);

figure(2)
vl_plotframe(f2);


[matches ,scores] = vl_ubcmatch(d1,d2);


x1 = [f1(1,matches(1,:));f1(2,matches(1,:))]; 
x2 = [f2(1,matches(2,:));f2(2,matches(2,:))];


perm = randperm(size(matches ,2)); figure;
imagesc([im1 im2]);
hold on;
plot([x1(1,perm(1:10)); x2(1,perm(1:10))+size(im1,2)], ... 
    [x1(2,perm(1:10)); x2(2,perm(1:10))],'-');
hold off;

save('for_comp_e5.mat', 'x1', 'x2');

%%% Answer
% How many of the matches appear to be correct?
% -> All of them actually.
