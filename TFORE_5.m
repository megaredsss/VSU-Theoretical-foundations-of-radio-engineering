clear all; close all; clc;
w=2*pi*10^4;
disp(w);
z1 = j*w*1.2*10^(-3); z2 = 9.5*10^3; z3 = 1/(j*w*2*10^(-9)); z4 = j*w*2.2*10^(-3); z5 = 12*10^3; 
z11 = z1+z2; z12 = -z2; z13 = 0; z21 = -z2; z22 = z2+z3+z4; z23 = -z4; z31 = 0; z32 = -z4; z33 = z4+z5;
Z = [z11 z12 z13;z21 z22 z23;z31 z32 z33];

e1 = 1.15*exp(j*pi/5);e2 = 2.2*exp(j*(-1)*pi/4);e3 = 0.9*exp(j*(-1)*pi/3);
E1 = e1-e2; E2 = e2+e3; E3=e3;
U = [E1; E2; E3;]

i1 = 5; i3 = 6; i5 = 7;
I = [i1; i3; i5;];
disp(Z^-1);
Z1=Z^-1;
Z2=Z1*U;
disp(Z);
disp('  ');
disp(U);
disp('  ');
disp(Z2);
