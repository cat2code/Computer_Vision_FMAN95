clear;clc;


e = 1.602*1e-19;
me = 9.11*1e-31;    
E0 = 231.94*1e9;    %%% kan vara fel
lambda = 800*1e-9;  % 800 nm
flen = 50*1e-2;     % 50 cm
D = 8*1e-3;         % 8 mm

d = 4*lambda*flen/pi/D;

P = 5*1e-3/22/1e-15; %peak power


A = pi*(d/2)^2;

Ipeak = P/A;


c = 299792458;
f = c/800/1e-9

w0 = f*2*pi;   % borde verkligen vara rätt

Up = (e^2*E0^2)/(4*me*w0^2)

