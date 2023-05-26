% Robot Controls - Robust control law for RRbot Manipulator
% Author: Krutarth Trivedi | ktrivedi@wpi.edu

clear; 
clc; 
close all;

% physical parameters of the robot
m1 = 1; r1 = 0.45; l1 = 1; I1 = 0.084;
m2 = 1; r2 = 0.45; l2 = 1; I2 = 0.084;
g = 9.81;

% define symbols
syms theta1 theta2 theta1_dot theta2_dot theta1_ddot theta2_ddot 'real'

%--- Initial Setup - (a) - Cubic Trajectory Generation ------%
syms t;
timeVec = [1; t; t^2; t^3];
timeVec_dot = diff(timeVec);
timeVec_ddot = diff(timeVec_dot);

t0 = 0; 
tf = 10;
theta1_t0 = pi;
theta1_tf = 0;
theta2_t0 = pi/2;
theta2_tf = 0;
theta1_dot_t0 = 0;
theta1_dot_tf = 0;
theta2_dot_t0 = 0;
theta2_dot_tf = 0;

timeMat = [1, t0, t0^2, t0^3;
           0, 1, 2*t0, 3*t0^2;
           1, tf, tf^2, tf^3;
           0, 1, 2*tf, 3*tf^2];

theta1ConfigVec = [theta1_t0; theta1_dot_t0; theta1_tf; theta1_dot_tf];
theta2ConfigVec = [theta2_t0; theta2_dot_t0; theta2_tf; theta2_dot_tf];

theta1CoefficientVec = pinv(timeMat)*theta1ConfigVec;
theta2CoefficientVec = pinv(timeMat)*theta2ConfigVec;

q1 = theta1CoefficientVec' * timeVec;
q2 = theta2CoefficientVec' * timeVec;

q1_dot = theta1CoefficientVec' * timeVec_dot;
q2_dot = theta2CoefficientVec' * timeVec_dot;

q1_ddot = theta1CoefficientVec' * timeVec_ddot;
q2_ddot = theta2CoefficientVec' * timeVec_ddot;

q_desired = [q1;q2];
qdot_desired = [q1_dot;q2_dot];
qddot_desired = [q1_ddot;q2_ddot];
Q_desired = [q_desired; qdot_desired];

%--------- Inital Setup - (b) - Manipulator Equation Form ---------%
% Nominal values as given in the problem statement
m1_hat = 0.75; m2_hat = 0.75;
I1_hat = 0.063; I2_hat = 0.063;

a = I1_hat + I2_hat + m1_hat*r1^2 + m2_hat*(l1^2 + r2^2);
b = m2_hat*l1*r2;
d = I2_hat + m2_hat*r2^2;

Mmat= [a+2*b* cos(theta2), d+b* cos(theta2); d+b* cos(theta2), d];
Cmat= [-b* sin(theta2)*theta2_dot, -b* sin(theta2)*(theta1_dot+theta2_dot); b* sin(theta2)*theta1_dot,0];
Gmat= [-m1_hat*g*r1* sin(theta1)-m2_hat*g*(l1* sin(theta1)+r2* sin(theta1+theta2)); -m2_hat*g*r2* sin(theta1+theta2)];

% ------- Robust Inverse Dynamics Control --------%
A = [0,0,1,0;
    0,0,0,1;
    0,0,0,0;
    0,0,0,0];

B = [0,0;
    0,0;
    1,0;
    0,1];

desiredEigenValues = [-1,-1,-2,-2];
K = place(A, B, desiredEigenValues);

Acl = A - B*K;
Q = eye(4);     
P = lyap(Acl',Q);
phi = 0.01;         %boundary layer - tunable Parameter
rho = 10;           %upper bound - tunable Parameter

% -------- Simulation and Plotting ---------%

% simulate the system for 10 sec for given control inputs using ODE45
T = 10;
y0 = [deg2rad(200), deg2rad(125), 0 ,0];
[t,y] = ode45(@ode_rrbot, 0:0.1:T, y0);

%reconstruct the desired trajectories and control input
q1_reconstruct = [];
q2_reconstruct = [];
q1_dot_reconstruct = [];
q2_dot_reconstruct= [];
q1_ddot_reconstruct = [];
q2_ddot_reconstruct= [];
t1 = [];
t2 = [];

for i = 1:size(t)
    q1_reconstruct(i) = double(subs(q1, t(i,:)));
    q2_reconstruct(i) = double(subs(q2, t(i,:)));
    q1_dot_reconstruct(i) = double(subs(q1_dot, t(i,:)));
    q2_dot_reconstruct(i) = double(subs(q2_dot, t(i,:)));
    q1_ddot_reconstruct(i) = double(subs(q1_ddot, t(i,:)));
    q2_ddot_reconstruct(i) = double(subs(q2_ddot, t(i,:)));

    error = [y(i,1);y(i,2);y(i,3);y(i,4)] - [q1_reconstruct(i);q2_reconstruct(i);q1_dot_reconstruct(i);q2_dot_reconstruct(i)];
    nr = error'*P*B;
    norm_nr = norm(nr);
    if phi > 0
        if phi >= norm_nr
            vr = -rho*(nr/phi);
        elseif phi < norm_nr
            vr = -rho*(nr/norm_nr);
        end
    else
        if norm_nr ~= 0
            vr = -rho*(nr/norm_nr);
        else
            vr = [0,0];
        end
    end

    V = -K*error + [q1_ddot_reconstruct(i);q2_ddot_reconstruct(i)] + vr';
    Mmat = double(subs(Mmat,y(i,2)));
    Cmat = double(subs(Cmat,[theta2,theta1_dot,theta2_dot],[y(i,2),y(i,3),y(i,4)]));
    Gmat = double(subs(Gmat,[theta1,theta2],[y(i,1),y(i,2)]));
    Tau = Mmat*V + Cmat*[y(i,3);y(i,4)] + Gmat;
    t1(i) = Tau(1);
    t2(i) = Tau(2);
end

% plot the trajectories
figure('Name','Trajectories', 'NumberTitle','off');
subplot(3,2,1)
plot(t,(y(:,1)),'b');
title('theta1')
xlabel('T');
ylabel('rad');
hold on;
plot(t,q1_reconstruct,'r');
legend({'Current','Desired'})

subplot(3,2,2)
plot(t,(y(:,2)),'b')
title('theta2')
xlabel('T');
ylabel('rad');
hold on;
plot(t,q2_reconstruct,'r');
legend({'Current','Desired'})

subplot(3,2,3)
plot(t,(y(:,3)),'b')
title('theta1-dot')
xlabel('T');
ylabel('rad/s');
hold on;
plot(t,q1_dot_reconstruct,'r');
legend({'Current','Desired'})

subplot(3,2,4);
plot(t,(y(:,4)),'b');
title('theta2-dot')
xlabel('T');
ylabel('rad/s');
hold on;
plot(t,q2_dot_reconstruct,'r');
legend({'Current','Desired'})

subplot(3,2,5);
plot(t,t1,'b');
title('t1')
xlabel('s');
ylabel('Nm');

subplot(3,2,6);
plot(t,t2,'b');
title('t2')
xlabel('s');
ylabel('Nm');

fprintf('Eigenvalues are selected such that: \n\n');
fprintf('Torque1: %.3f < u1 < %.3f \n\n',min(t1),max(t1));
fprintf('Torque2: %.3f < u2 < %.3f \n\n', min(t2),max(t2));