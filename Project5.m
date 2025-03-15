clear; close all; clc;

syms l1 l2 t1 t2 t3 m1 m2 dt1 dt1sq dt2 dt2sq ddt1 ddt2 g

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
%% Newton-Euler formulation
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
tau1 = simplify(transpose(n1) * [0; 0; 1]); % 
tau2 = simplify(transpose(n2) * [0; 0; 1]); % 
Tau = [tau1; tau2]; % Final Equation of motion at certain joint velocity and accelration, what will be joint torqe and acceleration

%% Lagrange formulation
%%







%% Design the Trajectory in the Joint Space

% Design the ViaPoints
initial = [0.1, 0];
final = [0.9, 0];
ViaPoints = [initial; final];

% RR1.plot([.1 0 0], 'workspace', [-1 1 -1 1 0 2]) %FIXME: Should it be three values in RR1

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
    
    % These calc angle, velocity, accel, jerk
    [t1q,t1qd,t1qdd,t1qddd] = linearParabolicBlendTrajectory(t0,tf,q0(1),qf(1),0,5,N);
    [t2q,t2qd,t2qdd,t2qddd] = linearParabolicBlendTrajectory(t0,tf,q0(2),qf(2),0,5,N);

    % Store results in the results matrix (concatenate col of 100)
    LPB_results(n:n+N-1, :) = [t1q(:), t2q(:), t1qd(:), t2qd(:), t1qdd(:), t2qdd(:), t1qddd(:), t2qddd(:)];
    n = n+N;
    disp("LPB Joint")
end

% %% Animate Joint Space Trajectory

%     % Must redo the config of the Links
%     L(1) = Link('revolute','d', 0, 'a', 0, 'alpha', 0 ,'modified');
%     L(2) = Link('revolute','d', 0, 'a', L1, 'alpha', 0 ,'modified');
%     L(3) = Link('revolute','d', 0, 'a', L2, 'alpha', 0 ,'modified'); % Needed for modified DH Table
%     RR1 = SerialLink(L, 'name', '2D-1-RR');

%     JointTrajectory = LPB_results;
%     t1 = JointTrajectory(:,1); % what joint angles you want to achieve
%     t2 = JointTrajectory(:,2);

%     % % Initialize VideoWriter
%     % JS = VideoWriter('JS_Traj_Finally.mp4', 'MPEG-4');
%     % open(JS);

%     figure()

%     for i = 1:height(JointTrajectory)
%         % Plot robot configuration for current time step
%         RR1.plot([t1(i) t2(i) 0], 'workspace', [-1 1 -1 1 0 2]) %FIXME: Should it be three values in RR1
%         xlabel('X'); ylabel('Y'); zlabel('Z');
%         title(['Joint Space: Frame ' num2str(i)]);

%         % % Capture current frame
%         % frame = getframe(gcf);
        
%         % % Write current frame to video
%         % writeVideo(JS, frame);
        
%         % Display progress
%         disp(['Running frame ' num2str(i)]);
%     end

%     % % Close VideoWriter
%     % close(JS);
%     % disp('Animation saved as joint_space_trajectory_animation.mp4');

% % End of Animation
%% Calculation with Gravity

for n = 1:length(t1q)
    t1 = t1q(n);
    dt1 = t1qd(n);
    ddt1 = t1qdd(n);
    t2 = t2q(n);
    dt2 = t1qd(n);
    ddt2 = t1qdd(n);

    [Tau, tauInertial, tauCent, tauCor, tauGrav] = calcJointTorque(L1, t1, dt1, ddt1, L2, t2, dt2, ddt2, 9.81);

    % Store the joint torques in TauResults
    TauResults(n, :) = [Tau(1), Tau(2), tauInertial(1), tauInertial(2), tauCent(1), tauCent(2), tauCor(1), tauCor(2), tauGrav(1), tauGrav(2)];
end

%% Calculation without Gravity
for n = 1:length(t1q)
    t1 = t1q(n);
    dt1 = t1qd(n);
    ddt1 = t1qdd(n);
    t2 = t2q(n);
    dt2 = t1qd(n);
    ddt2 = t1qdd(n);

    [Tau, tauInertial, tauCent, tauCor, tauGrav] = calcJointTorque(L1, t1, dt1, ddt1, L2, t2, dt2, ddt2, 0);

    % Store the joint torques in TauResults
    TauResultsNoGrav(n, :) = [Tau(1), Tau(2), tauInertial(1), tauInertial(2), tauCent(1), tauCent(2), tauCor(1), tauCor(2), tauGrav(1), tauGrav(2)];
end

%% Graph
% Define time vector (10 points over 4 seconds)
t = linspace(0, 4, 10);

% % Theta 1 Graphs (with Gravity)
%     figure;
%     % Create a 3x1 subplot layout
%     subplot(3, 1, 1);
%     plot(t, t1q, '-o'); % Plot data1 in the first subplot
%     xlabel('Time (s)');
%     ylabel('Angle (rad)');

%     subplot(3, 1, 2);
%     plot(t, t1qd, '-s'); % Plot data2 in the second subplot
%     xlabel('Time (s)');
%     ylabel('Velocity (rad/s)');

%     subplot(3, 1, 3);
%     plot(t, t1qdd, '-^'); % Plot data3 in the third subplot
%     xlabel('Time (s)');
%     ylabel('Accel. (rad/s^2)');

%     % Adjust subplot spacing
%     sgtitle('Joint Angles'); % Overall title for all subplots
% % End Theta 1 Graphs

% Theta Graphs (with Gravity)
    figure;
    % Create a 3x1 subplot layout
    subplot(3, 1, 1);
    plot(t, t1q, 'b', 'LineWidth', 2, 'DisplayName', 'theta1'); % Plot data1 in the first subplot
    hold on;
    plot(t, t2q, 'r', 'LineWidth', 2, 'DisplayName', 'theta2'); % Plot data1 in the first subplot
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('show', 'Location', 'best');
    title('Joint Angles');


    subplot(3, 1, 2);
    plot(t, t1qd, 'b', 'LineWidth', 2, 'DisplayName', 'omega1'); % Plot data2 in the second subplot
    hold on;
    plot(t, t2qd, 'r', 'LineWidth', 2, 'DisplayName', 'omega2'); % Plot data2 in the second subplot
    xlabel('Time (s)');
    ylabel('Velocity (rad/s)');
    legend('show', 'Location', 'best');
    title('Joint Velocity');


    subplot(3, 1, 3);
    plot(t, t2qdd, 'b', 'LineWidth', 2, 'DisplayName', 'sigma1'); % Plot data3 in the third subplot
    hold on;
    plot(t, t2qdd, 'r', 'LineWidth', 2, 'DisplayName', 'sigma2'); % Plot data3 in the third subplot
    xlabel('Time (s)');
    ylabel('Accel. (rad/s^2)');
    legend('show', 'Location', 'best');
    title('Joint Acceleration');


% End of Theta Graphs

% Torque Graphs (with Gravity)
    figure;

    % Create a 2x1 subplot layout
    subplot(2, 1, 1);
    plot(t, TauResults(:,3), '-b', 'LineWidth', 2,'DisplayName', 'Inertial');
    hold on;
    plot(t, TauResults(:,5), '-r', 'LineWidth', 2,'DisplayName', 'Centrifugal');
    plot(t, TauResults(:,7), '-y', 'LineWidth', 2,'DisplayName', 'Corolis');
    plot(t, TauResults(:,9), '-m', 'LineWidth', 2,'DisplayName', 'Gravity');
    plot(t, TauResults(:,1), '-g', 'LineWidth', 2,'DisplayName', 'Total');

    xlabel('Time (s)');
    ylabel('Torque (Nm)');
    title('Joint 1 Torques (with Gravity)');
    grid on;
    legend('show', 'Location', 'best');


    subplot(2, 1, 2);
    plot(t, TauResults(:,4), '-b', 'LineWidth', 2,'DisplayName', 'Inertial');
    hold on;
    plot(t, TauResults(:,6), '-r', 'LineWidth', 2,'DisplayName', 'Centrifugal');
    plot(t, TauResults(:,8), '-y', 'LineWidth', 2,'DisplayName', 'Corolis');
    plot(t, TauResults(:,10), 'm', 'LineWidth', 2,'DisplayName', 'Gravity');
    plot(t, TauResults(:,2), '-g', 'LineWidth', 2,'DisplayName', 'Total');
    xlabel('Time (s)');
    ylabel('Torque (Nm)');
    title('Joint 2 Torques (with Gravity)');
    grid on;
    legend('show');

% End of Torque Graphs

% Theta Graphs (WITHOUT Gravity)
    figure;
    % Create a 3x1 subplot layout
    subplot(3, 1, 1);
    plot(t, t1q, 'b', 'LineWidth', 2, 'DisplayName', 'theta1'); % Plot data1 in the first subplot
    hold on;
    plot(t, t2q, 'r', 'LineWidth', 2, 'DisplayName', 'theta2'); % Plot data1 in the first subplot
    xlabel('Time (s)');
    ylabel('Angle (rad)');
    legend('show', 'Location', 'best');
    title('Joint Angles');


    subplot(3, 1, 2);
    plot(t, t1qd, 'b', 'LineWidth', 2, 'DisplayName', 'omega1'); % Plot data2 in the second subplot
    hold on;
    plot(t, t2qd, 'r', 'LineWidth', 2, 'DisplayName', 'omega2'); % Plot data2 in the second subplot
    xlabel('Time (s)');
    ylabel('Velocity (rad/s)');
    legend('show', 'Location', 'best');
    title('Joint Velocity');


    subplot(3, 1, 3);
    plot(t, t2qdd, 'b', 'LineWidth', 2, 'DisplayName', 'sigma1'); % Plot data3 in the third subplot
    hold on;
    plot(t, t2qdd, 'r', 'LineWidth', 2, 'DisplayName', 'sigma2'); % Plot data3 in the third subplot
    xlabel('Time (s)');
    ylabel('Accel. (rad/s^2)');
    legend('show', 'Location', 'best');
    title('Joint Acceleration');

% End of Theta Graphs

% Torque Graphs (WITHOUT Gravity)
    figure;

    % Create a 2x1 subplot layout
    subplot(2, 1, 1);
    plot(t, TauResultsNoGrav(:,3), '-b', 'LineWidth', 2,'DisplayName', 'Inertial');
    hold on;
    plot(t, TauResultsNoGrav(:,5), '-r', 'LineWidth', 2,'DisplayName', 'Centrifugal');
    plot(t, TauResultsNoGrav(:,7), '-y', 'LineWidth', 2,'DisplayName', 'Corolis');
    plot(t, TauResultsNoGrav(:,9), '-m', 'LineWidth', 2,'DisplayName', 'Gravity');
    plot(t, TauResultsNoGrav(:,1), '-g', 'LineWidth', 2,'DisplayName', 'Total');
    xlabel('Time (s)');
    ylabel('Torque (Nm)');
    title('Joint 1 Torques (without Gravity)');
    grid on;
    legend('show', 'Location', 'best');

    subplot(2, 1, 2);
    plot(t, TauResultsNoGrav(:,4), '-b', 'LineWidth', 2,'DisplayName', 'Inertial');
    hold on;
    plot(t, TauResultsNoGrav(:,6), '-r', 'LineWidth', 2,'DisplayName', 'Centrifugal');
    plot(t, TauResultsNoGrav(:,8), '-y', 'LineWidth', 2,'DisplayName', 'Corolis');
    plot(t, TauResultsNoGrav(:,10), 'm', 'LineWidth', 2,'DisplayName', 'Gravity');
    plot(t, TauResultsNoGrav(:,2), '-g', 'LineWidth', 2,'DisplayName', 'Total');
    xlabel('Time (s)');
    ylabel('Torque (Nm)');
    title('Joint 2 Torques (without Gravity)');
    grid on;
    legend('show');

% End of Torque Graphs

% Torque Comparison
    figure;

    % Create a 2x1 subplot layout
    plot(t, TauResultsNoGrav(:,1)./TauResults(:,1), '-b', 'LineWidth', 2,'DisplayName', 'Joint 1');
    hold on;
    plot(t, TauResultsNoGrav(:,2)./TauResults(:,2), '-r', 'LineWidth', 2,'DisplayName', 'Joint 2');
    xlabel('Time (s)');
    ylabel('Torque (Nm)');
    title('Torque Comparison');
    grid on;
    legend('show', 'Location', 'best');
% End of Torque Comparison

%% OTHER

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


% Tau is torque from end effector
% N1 - from joint
% n1 - from com

% N1 is dynamic
% n1 combine static and dynamic
% Tau - just the third row just in that direction, other rows are in the mechanism structure


% +N1 is net torqye from statis and dynamic, exert on EE and due to motion
% n1 - just static
% Tau1


% TODO: Create a function to make the whole thing again while inputting values from the trajectory design
% TODO: Redo into the fact of Inertia Values

%% Problem 2
rho = 2710;

d = 0.4;
vA = 0.1^3 - (0.01^2 * pi * 0.1); vB = (0.05/2)^2 * pi * 0.7; vC = vA;
mA = vA * rho; mB = vB * rho; mC = vC * rho;
% body A
I_A = mA/12*[0.1^2+0.1^2 0 0; 0 0.1^2+0.1^2 0; 0 0 0.1^2+0.1^2] - [(mA/12)*(3*0.01^2+0.1^2) 0 0; 0 (mA/12)*(3*0.01^2+0.1^2) 0; 0 0 (mA*0.01^2)/2];
I_Acm = I_A+mA*[0 0 0; 0 d^2 0; 0 0 d^2];
% body B
I_Bcm = [(mB*0.025^2)/2 0 0; (mB/12)*(3*0.025^2+0.8^2) 0 0; 0 (mB/12)*(3*0.025^2+0.8^2) 0];
% body C
I_C = mC/12*[0.1^2+0.1^2 0 0; 0 0.1^2+0.1^2 0; 0 0 0.1^2+0.1^2] - [(mC/12)*(3*0.01^2+0.1^2) 0 0; 0 (mC/12)*(3*0.01^2+0.1^2) 0; 0 0 (mC*0.01^2)/2];
R = rotx(-45);
I_C = R*I_C*transpose(R);
I_Ccm = I_C+mC*[0 0 0; 0 d^2 0; 0 0 d^2];
I_tot = I_Acm + I_Bcm + I_Ccm;
fprintf('Tensor of Inertia:\n');
disp(I_tot);