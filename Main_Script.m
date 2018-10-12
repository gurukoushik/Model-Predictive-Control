%Setting up data required for implementing the MPC, state estimation and
%plant models.
clear all
load('linssmodel.mat');
load('DATA_REQD.mat');
fcc_parameters
C=C_new;
Hc=4; %Control Horizon
Hp=8; %Prediction Horizon
weights=[10 10 1000]; %weights for phi matrix in the objective function
yref=[0.03;0.25;375];%Set point
y_traj=yref; %Reference trajectory
Q=1*eye(size(A,1)); %Process covariance
R=0.09*eye(size(C,1)); %Measurement covariance
P0=100*eye(size(A,1)); %Estimation covariance initial
Xinit=X0; %initial state estimate
Pinit=P0;
YY=C*Xinit; %initial measurement estimate using state estimate 
for i=1:100
%y_traj=yref*(1-exp(-5*i+5)); %Exponential reference trajecory
%y_traj=yref*i/100; %Straight line reference trajectory
inputFromMPC=objectivemin(A,B,C,Hc,Hp,Xinit,y_traj,weights); %Calling function to return optimal inputs
Ynew=Plant(i,YY,Time,inputFromMPC); %Passing optimal inputs to plant model and getting new measurements of state 
Xup=A*Xinit+B*inputFromMPC; %Propogation step for Kalman filter
Pup=A*Pinit*transpose(A)+Q; %Covariance propogation for Kalman filter
K=Pup*transpose(C)/(C*Pup*transpose(C)+R); %Calculating Kalman Gain
Xcor=Xup+K*(transpose(Ynew)-C*Xup); %Correction step for Kalman filter
Pcor=Pup-K*C*Pup; %Correction step for covariance in Kalman filter
Xest(:,i)=Xcor;
Ymeas(:,i)=C*Xest(:,i);
Xinit=Xcor;
Pinit=Pcor;
YY=Ynew;
inputMPCone(i)=inputFromMPC(1);
inputMPCtwo(i)=inputFromMPC(2);
end
%Plot Y and U
figure(1);
plot(Ymeas(1,:));
title('Concentration 1 VS Time');
xlabel('Time');
ylabel('Concentration 1')
figure(2);
plot(Ymeas(2,:));
title('Concentration 2 VS Time');
xlabel('Time');
ylabel('Concentration 2')
figure(3);
plot(Ymeas(3,:));
title('Temperature VS Time');
xlabel('Time');
ylabel('Temperature')
figure(4);
stairs(inputMPCone);
title('dFa VS Time');
xlabel('Time');
ylabel('dFa')
figure(5);
stairs(inputMPCtwo);
title('dFsc VS Time');
xlabel('Time');
ylabel('dFsc')