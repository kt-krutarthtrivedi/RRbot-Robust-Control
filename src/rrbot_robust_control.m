% Robot Controls - Robust control law for RRbot Manipulator
% Author: Krutarth Trivedi | ktrivedi@wpi.edu

clear; close; clc;

% ROS Setup
rosinit;

% physical parameters of the robot
m1 = 1; r1 = 0.45; l1 = 1; I1 = 0.084;
m2 = 1; r2 = 0.45; l2 = 1; I2 = 0.084;
g = 9.81;

j1_effort = rospublisher('/rrbot/joint1_effort_controller/command');
j2_effort = rospublisher('/rrbot/joint2_effort_controller/command');
JointStates = rossubscriber('/rrbot/joint_states');

tau1 = rosmessage(j1_effort);
tau2 = rosmessage(j2_effort);

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

client = rossvcclient('/gazebo/set_model_configuration');
req = rosmessage(client);
req.ModelName = 'rrbot';
req.UrdfParamName = 'robot_description';
req.JointNames = {'joint1','joint2'};
req.JointPositions = [deg2rad(200), deg2rad(125)];
resp = call(client,req,'Timeout',3);

tic;
t = 0;

sample_theta1 = [];
sample_theta2 = [];
sample_theta1_dot = [];
sample_theta2_dot = [];
sample_time = [];
sample_tau1 = [];
sample_tau2 = [];
sample_q1 = [];
sample_q2 = [];
sample_q1_dot = [];
sample_q2_dot = [];
sample_q1_ddot = [];
sample_q2_ddot = [];

m1_hat = 0.75; m2_hat = 0.75;
I1_hat = 0.063; I2_hat = 0.063;

% ------- Robust Inverse Dynamics Control --------%
B = [0,0;
    0,0;
    1,0;
    0,1];

K = [2.0000,0, 3.0000, 0;
     0, 2.0000, 0, 3.0000];

P = [1.2500,        0,    0.2500,         0;
         0,    1.2500,         0,    0.2500;
    0.2500,         0,    0.2500,         0;
         0,    0.2500,         0,    0.2500];

phi = 0.1;         %boundary layer - tunable Parameter
rho = 10;          %upper bound - tunable Parameter

while(t < 10)

    t = toc;
    % read the joint states
    jointData = receive(JointStates);

    % implementing the trajectory feedback controller  
    theta1 = jointData.Position(1,1);
    theta1 = wrapTo2Pi(theta1);
    theta2 = jointData.Position(2,1);
    theta2 = wrapTo2Pi(theta2);
    theta1_dot = jointData.Velocity(1,1);
    theta2_dot = jointData.Velocity(2,1);

    jointValues = [theta1; theta2; theta1_dot; theta2_dot];

    %--- Generate cubic polynomial trajectories for both the joints ------%
     
    q1 = (pi*t^3)/500 - (3*pi*t^2)/100 - (6189958033024885*t)/10141204801825835211973625643008 + pi;
    q2 = (pi*t^3)/1000 - (3*pi*t^2)/200 - (6189958033024885*t)/20282409603651670423947251286016 + pi/2;
    q1_dot = (3*pi*t^2)/500 - (3*pi*t)/50 - 6189958033024885/10141204801825835211973625643008;
    q2_dot = (3*pi*t^2)/1000 - (3*pi*t)/100 - 6189958033024885/20282409603651670423947251286016;
    q1_ddot = (3*pi*t)/250 - (3*pi)/50; 
    q2_ddot = (3*pi*t)/500 - (3*pi)/100;

    q_desired = [q1;q2];
    qdot_desired = [q1_dot;q2_dot];
    qddot_desired = [q1_ddot;q2_ddot];
    Q_desired = [q_desired; qdot_desired];
    
    %--------- Manipulator Equation Form ---------%    
    a = I1_hat + I2_hat + m1_hat*r1^2 + m2_hat*(l1^2 + r2^2);
    b = m2_hat*l1*r2;
    d = I2_hat + m2_hat*r2^2;
    
    Mmat= [a+2*b* cos(theta2), d+b* cos(theta2); d+b* cos(theta2), d];
    Cmat= [-b* sin(theta2)*theta2_dot, -b* sin(theta2)*(theta1_dot+theta2_dot); b* sin(theta2)*theta1_dot,0];
    Gmat= [-m1_hat*g*r1* sin(theta1)-m2_hat*g*(l1* sin(theta1)+r2* sin(theta1+theta2)); -m2_hat*g*r2* sin(theta1+theta2)];
    
    error = jointValues - Q_desired;
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
    
    %--------------- Virtual control input ----------------------%
    V = -K*error + qddot_desired + vr';
    
    %-------------- The overall control law ---------------%
    Tau = Mmat*V + Cmat*[theta1_dot; theta2_dot] + Gmat;
 
    tau1.Data = Tau(1);
    tau2.Data = Tau(2);
    
    send(j1_effort,tau1);
    send(j2_effort,tau2);

    % sample the time, joint state values, and calculated torques here to be plotted at the end
    sample_time = [sample_time, t];
    sample_theta1 = [sample_theta1, theta1];
    sample_theta2 = [sample_theta2, theta2];
    sample_theta1_dot = [sample_theta1_dot, theta1_dot];
    sample_theta2_dot = [sample_theta2_dot, theta2_dot];
    sample_tau1 = [sample_tau1, Tau(1)];
    sample_tau2 = [sample_tau2, Tau(2)];
    sample_q1 = [sample_q1, q1];
    sample_q2 = [sample_q2, q2];
    sample_q1_dot = [sample_q1_dot, q1_dot];
    sample_q2_dot = [sample_q2_dot, q2_dot];
    sample_q1_ddot = [sample_q1_ddot, q1_ddot];
    sample_q2_ddot = [sample_q2_ddot, q2_ddot];
end

tau1.Data = 0;
tau2.Data = 0;

send(j1_effort,tau1);
send(j2_effort,tau2);

% disconnect from roscore
rosshutdown;

% plot the trajectories
figure('Name','State Trajectories', 'NumberTitle','off');
subplot(3,2,1)
plot(sample_time,sample_theta1,'b');
title('theta1')
xlabel('T');
ylabel('rad');
hold on;
plot(sample_time,sample_q1,'r');
legend({'Current','Desired'})

subplot(3,2,2)
plot(sample_time,sample_theta2,'b')
title('theta2')
xlabel('T');
ylabel('rad');
hold on;
plot(sample_time,sample_q2,'r');
legend({'Current','Desired'})

subplot(3,2,3)
plot(sample_time,sample_theta1_dot,'b')
title('theta1-dot')
xlabel('T');
ylabel('rad/s');
hold on;
plot(sample_time,sample_q1_dot,'r');
legend({'Current','Desired'})

subplot(3,2,4);
plot(sample_time,sample_theta2_dot,'b');
title('theta2-dot')
xlabel('T');
ylabel('rad/s');
hold on;
plot(sample_time,sample_q2_dot,'r');
legend({'Current','Desired'})

subplot(3,2,5);
plot(sample_time,sample_tau1,'b');
title('tau1')
xlabel('s');
ylabel('Nm');

subplot(3,2,6);
plot(sample_time,sample_tau2,'b');
title('tau2')
xlabel('s');
ylabel('Nm');