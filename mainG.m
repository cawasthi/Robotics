%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code performs trajectory optimization on a planar cart pole problem
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



% % Direct transcription
clear;
clc;


global nGrid dt

% Defining grid
nGrid = 80;
ti = 0;
tf = 2;
t(1,:) = linspace(ti,tf,nGrid);
dt = (t(end)-t(1))/(nGrid-1);
d = 1;
dmax = 2;
umax = 20;


% Guess to initialize the NLP
% x is state, q is position
tGuess = [0,tf];                % time corresponding to xGuess and uGuess
xGuess = [0,0,0,0;d,pi,0,0];    % initial and final value of state traj
uGuess = [0;0];                 % inital and final value of control traj


% Interpolate guess at the grid points for transcription
u_init(:,1)     = interp1(tGuess,uGuess,t,'linear');   % linear interpolation of control vector, initial jspace state traj
state_init(:,1:4) = interp1(tGuess,xGuess,t,'linear');   % linear interpolation of state vector, initial jspace control traj


% simple linear interpolation trajectory
uu = reshape(u_init',[nGrid,1]);
qq = reshape(state_init',[4*nGrid,1]);
x0 = [uu;qq];


A   = [];
b   = [];
Aeq = [];
beq = [];


lb = [];
ub = [];



%Options for fmincon
% options = optimoptions(@fmincon, 'TolFun', 1.0e-06,...
%                        'OptimalityTolerance',1.0e-06,'MaxIter', 10000,'Display','iter','MaxFunEvals', 150000, ...
%                        'StepTolerance',1.0e-06,'DiffMinChange', 0.001, 'Algorithm', 'sqp');

% options = optimoptions(@fmincon, 'SpecifyObjectiveGradient',true,'TolFun', 1.0e-06,'checkgradients',true,'FiniteDifferenceType','central',...
%                       'OptimalityTolerance',1.0e-06,'MaxIter', 10000,'Display','iter','MaxFunEvals', 200000, ...
%                       'StepTolerance',1.0e-06,'DiffMinChange', 0.001, 'Algorithm', 'sqp');
                   
      
options = optimoptions(@fmincon, 'SpecifyObjectiveGradient',true,'TolFun', 1.0e-06,...
                      'OptimalityTolerance',1.0e-06,'MaxIter', 10000,'Display','iter','MaxFunEvals', 200000, ...
                      'StepTolerance',1.0e-06,'DiffMinChange', 0.001, 'Algorithm', 'sqp');

% options = optimoptions(@fmincon, 'TolFun', 1.0e-06,'FiniteDifferenceType','forward',...
%                       'OptimalityTolerance',1.0e-06,'MaxIter', 10000,'Display','iter','MaxFunEvals', 100000, ...
%                       'StepTolerance',1.0e-06,'DiffMinChange', 0.001, 'Algorithm', 'sqp');
                   
% Solve the constrained optimization problem
optimal = fmincon(@pathObjWithGrad, x0, A, b, Aeq, beq, lb, ub, ...
              @dynConsG, options);


% optimal and guess linear position vs. time          
% subplot(3,1,1)
% plot(t,x0(nGrid+1:4:end-3),t,optimal(nGrid+1:4:end-3),'b.','MarkerSize',12)          
% ylim([0 1.5])
% 
% % optimal and guess angular position vs. time
% subplot(3,1,2)
% plot(t,x0(nGrid+2:4:end-2),t,optimal(nGrid+2:4:end-2),'b.','MarkerSize',12)
% ylim([-2 4])
% 
% % optimal and guess control vs. time
% subplot(3,1,3)
% plot(t,x0(1:nGrid),t, optimal(1:nGrid),'b.','MarkerSize',12)
% ylim([-20 11])

          
% % Optimal state and control trajectories          
q_opt = [optimal(nGrid+1:4:end-3),optimal(nGrid+2:4:end-2)];
qdot_opt = [optimal(nGrid+3:4:end-1),optimal(nGrid+4:4:end)];
u_opt = optimal(1:nGrid);


save('optimVarsCP.mat','q_opt','qdot_opt','u_opt','dt','nGrid')
          
          
% %%% Plots of joint positions, joint velocities and joint torques          
%           
% subplot(3,1,1)
% plot(t,optimal(3*nGrid+1:6:end-5),t,optimal(3*nGrid+2:6:end-4),t,optimal(3*nGrid+3:6:end-3))
% legend('q1','q2','q3')
% xlabel('t (in sec)')
% ylabel('joint position (in rad)')
% grid on
% 
% 
% subplot(3,1,2)
% plot(t,optimal(3*nGrid+4:6:end-2),t,optimal(3*nGrid+5:6:end-1),t,optimal(3*nGrid+6:6:end))
% legend('q1 dot','q2 dot','q3 dot')
% xlabel('t (in sec)')
% ylabel('joint velocity (in rad/s)')
% grid on
%  
%           
% subplot(3,1,3)
% plot(t,optimal(1:3:3*nGrid),t,optimal(2:3:3*nGrid-1),t,optimal(3:3:3*nGrid)-2)
% legend('u1','u2','u3')
% xlabel('t (in sec)')
% ylabel('joint torque (in N-m)')
% grid on





