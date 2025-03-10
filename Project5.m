% clear; close all; clc;

syms l1 l2 t1 t2 t3 m1 m2 dt1 dt2 ddt1 ddt2 g

L(1) = Link('revolute','d', 0, 'a', 0, 'alpha', 0 ,'modified');
L(2) = Link('revolute','d', 0, 'a', l1, 'alpha', 0 ,'modified');
L(3) = Link('revolute','d', 0, 'a', l2, 'alpha', 0 ,'modified'); % Needed for modified DH Table

RR1 = SerialLink(L, 'name', '2D-1-RR');
th = [t1 t2 0];

T0_1 = RR1.A([1], th);
T1_2 = RR1.A([2], th);
T2_T = RR1.A([3], th);
T0_T = RR1.A([1 2 3], th);
%%
[R0_1, P0_1] = tr2rt(T0_1); R1_0 = transpose(R0_1);
[R1_2, P1_2] = tr2rt(T1_2); R2_1 = transpose(R1_2);
[R2_T, P2_T] = tr2rt(T2_T); RT_2 = transpose(R2_T);
[R0_T, P0_T] = tr2rt(T0_T);
%%
% Center of Mass Location (l1/2)
PC1 = [l1/2; 0; 0]; % Fill this out
PC2 = [l2/2; 0; 0]; % Fill this out

IC1 = (1/12) * m1 * l1^2 * [0 0 0; 0 1 0; 0 0 1]; % Moment of Inertia
IC2 = (1/12) * m2 * l2^2 * [0 0 0; 0 1 0; 0 0 1];

f3 = zeros(3,1); % External force at EE
n3 = zeros(3,1); % External Torque at EE

% values at base frame
w0 = zeros(3,1); % rot vel
wd0 = zeros(3,1); % rot accel

v0 = zeros(3,1); % linear vel
vd0 = [0; g; 0]; % linear accel (gravity effects)
%% outward iteration to compete velocities and accelerations
% z is always 001; ddt is thea double dot; dt is theta dot

% i = 0
w1 = R1_0 * w0 + dt1 * [0; 0; 1]; % Fill this out
wd1 = R1_0 * wd0 + R1_0 * cross(w0, dt1 * [0; 0; 1]) + ddt1 * [0; 0; 1]; % THIS IS WRONG PARANTEHSIS

vd1 = R1_0 * (cross(wd0,P0_1) + cross(w0,cross(w0,P0_1)) + vd0); % Fill this out
vcd1 = cross(wd1,PC1) + cross(w1, cross(w1,PC1)) + vd1; % Fill this out

F1 = m1 * vcd1; % Fill this out
N1 = IC1 * wd1 + cross(w1,IC1*w1); % Fill this out

% i = 1
w2 = R2_1 * w1 + dt2 * [0; 0; 1]; % Fill this out
wd2 = R2_1 * wd1 + R2_1 * cross(w1, dt2 * [0; 0; 1])+ ddt2 * [0; 0; 1;]; % Fill this out

vd2 = R2_1 * (cross(wd1,P1_2) + cross(w1,cross(w1,P1_2)) + vd1); % Fill this out
vcd2 = cross(wd2,PC2) + cross(w2, cross(w2,PC2)) + vd2; % Fill this out

F2 = m2 * vcd2; % Fill this out
N2 = IC2 * wd2 + cross(w2,IC2*w2); % Fill this out
%% inward iteration to compute forces and torques
% i = 2
f2 = R2_T * f3 + F2; % Fill this out % Uppercase is due to dynamics
n2 = N2 + R2_T * n3 + cross (PC2, F2) + cross (P2_T, R2_T * f3); % Fill this out
% i = 1
f1 = R1_2 * f2 + F1; % Fill this out
n1 = N1 + R1_2 * n2 + cross (PC2, F1) + cross (P1_2, R1_2 * f2); % Fill this out
%% Derive final result
tau1 = transpose(n1) * [0; 0; 1]; % 
tau2 = transpose(n2) * [0; 0; 1]; % 
Tau = [tau1; tau2]; % Final Equation of motion at certain joint velocity and accelration, what will be joint torqe and acceleration
% tau 1 will be the torque or force of the value that we get

%2.7 grams al 6061
% solve conjugate issue with putting fake values in the sysm variables
% can also do system preference for doing decimal



% add 10 N on the f3 term, 

% inverse from .1 and .9 and do joint angle, will get angular for these

% remove gravity and compare what is different from dynamic equation

% Direct reference of how to do this in the reference notes do 2.7 grams/cubic cm for the density

% collect function will get the coefficients for the values

% intuitive robotics (3 rounds), nova medical

% qdd = .8726