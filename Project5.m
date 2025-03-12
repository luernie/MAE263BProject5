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

%% Design the Trajectory in the Joint Space

% Design the ViaPoints
initial = [0.1, 0];
final = [0.9, 0];
ViaPoints = [initial; final];

RR1.plot([.1 0 0], 'workspace', [-1 1 -1 1 0 2]) %FIXME: Should it be three values in RR1


L1 = 0.5; % Arm Lengths %FIXME:Can I assume these lengths? 
L2 = 0.5;

% Calculate IK for each point
for i = 1:height(ViaPoints)
    % Perform inverse kinematics calculation
    sol = IK(ViaPoints(i,:), L1, L2, 0);
    
    % Store results in the results matrix
    JointWaypoints(i, :) = sol;
%   disp("run")
end

%% Linear Paraoblic Blend for Joint Space Trajectory %TODO: Issue should be here t1 changeing, t2 constant like clock
t0 = 0; % initial time
tf = 4; % final time
N = 10; % number of points

% Loop through LPB multiple times
n = 1;
for i = 1:height(ViaPoints)-1
    % Enter the LPB for X,Y,Z,EEorient
    q0 = JointWaypoints(i,:);
    qf = JointWaypoints(i+1,:);
    [t1q,t1qd,t1qdd,t1qddd] = linearParabolicBlendTrajectory_inclass(t0,tf,q0(1),qf(1),0,50,N);
    [t2q,t2qd,t2qdd,t2qddd] = linearParabolicBlendTrajectory_inclass(t0,tf,q0(2),qf(2),0,50,N);

    % Store results in the results matrix (concatenate col of 100)
    LPB_results(n:n+N-1, :) = [t1q(:), t2q(:), t1qd(:), t2qd(:), t1qdd(:), t2qdd(:), t1qddd(:), t2qddd(:)];
    n = n+N;
    disp("LPB Joint")
end

%% Animate Joint Space Trajectory

% Must redo the config of the Links
L(1) = Link('revolute','d', 0, 'a', 0, 'alpha', 0 ,'modified');
L(2) = Link('revolute','d', 0, 'a', L1, 'alpha', 0 ,'modified');
L(3) = Link('revolute','d', 0, 'a', L2, 'alpha', 0 ,'modified'); % Needed for modified DH Table
RR1 = SerialLink(L, 'name', '2D-1-RR');

JointTrajectory = LPB_results;
t1 = JointTrajectory(:,1); % what joint angles you want to achieve
t2 = JointTrajectory(:,2);

% % Initialize VideoWriter
% JS = VideoWriter('JS_Traj_Finally.mp4', 'MPEG-4');
% open(JS);

figure()

for i = 1:height(JointTrajectory)
    % Plot robot configuration for current time step
    RR1.plot([t1(i) t2(i) 0], 'workspace', [-1 1 -1 1 0 2]) %FIXME: Should it be three values in RR1
    xlabel('X'); ylabel('Y'); zlabel('Z');
    title(['Joint Space: Frame ' num2str(i)]);

    % % Capture current frame
    % frame = getframe(gcf);
    
    % % Write current frame to video
    % writeVideo(JS, frame);
    
    % Display progress
    disp(['Saving frame ' num2str(i)]);
end

% % Close VideoWriter
% close(JS);
% disp('Animation saved as joint_space_trajectory_animation.mp4');

%%
q = linearParabolicBlendTrajectory(t0, tf, q0, qf, tb, qdd, N)%


X = 0:0.1:1; % X-axis positions from 0 to 1 meter in 0.1 meter increments
Y = 0;

p = [1,1];

s = 0; % elbow up value

sol = IK(p, L1, L2, s); % spits out two angles


%% Linear Paraoblic Blend for Task Space
t0 = 0; % initial time
tf = 5; % final time
N = 10; % number of points

% Initialize the Matrix
LPB_Task_results = zeros(N, 16);

% Loop through LPB multiple times
n = 1;
for i = 1:height(TaskSequence)-1
    % Enter the LPB for X,Y,Z,EEorient
    q0 = TaskSequence(i,:);
    qf = TaskSequence(i+1,:);
    [Xq,Xqd,Xqdd,Xqddd] = linearParabolicBlendTrajectory_inclass(t0,tf,q0(1),qf(1),0,50,N);
    [Yq,Yqd,Yqdd,Yqddd] = linearParabolicBlendTrajectory_inclass(t0,tf,q0(2),qf(2),0,50,N);
    [Zq,Zqd,Zqdd,Zqddd] = linearParabolicBlendTrajectory_inclass(t0,tf,q0(3),qf(3),0,50,N);
    [EEq,EEqd,EEqdd,EEqddd] = linearParabolicBlendTrajectory_inclass(t0,tf,q0(4),qf(4),0,50,N);

    % Store results in the results matrix (concatenate col of 100)
    LPB_Task_results(n:n+N-1, :) = [Xq(:), Yq(:), Zq(:), EEq(:), Xqd(:), Yqd(:), Zqd(:), EEqd(:), Xqdd(:), Yqdd(:), Zqdd(:), EEqdd(:), Xqddd(:), Yqddd(:), Zqddd(:), EEqddd(:)];
    n = n+N;
    disp("LPB Task")
end






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