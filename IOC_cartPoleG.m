% % IOC for non-planar robot

% clear workspace
clear;
clc;
syms x1 x2 x3 x4 u1


% % Forward Control % %

% load optimal state and control trajectories from nonlinear controller

load('optimVarsCP.mat','q_opt','qdot_opt','u_opt','dt','nGrid')
dt = 2/(nGrid-1);

% load('q1_opt.mat','Xspline','Yspline')
% q1_opt = Yspline;
% time = Xspline';
% 
% load('q2_opt.mat','Xspline','Yspline')
% q2_opt = Yspline;
% 
% load('q1dot_opt.mat','Xspline','Yspline')
% q1dot_opt = Yspline;
% 
% load('q2dot_opt.mat','Xspline','Yspline')
% q2dot_opt = Yspline;
% 
% load('u_opt.mat','Xspline','Yspline')
% u_opt = Yspline;
% 
% 
% % %preprocessing
% q_opt = [q1_opt' q2_opt'];
% qdot_opt = [q1dot_opt' q2dot_opt'];
% u_opt = u_opt';
% 
% dt = 2/2302; 


% % Inverse control % %

% The Lagrangian of the original cost function is given bt c.'*phi
% c = cost vectors (needed to be found)
% phi = feature vector set (need to be guessed from domain knowledge)
%
% In our case, x = 4x1 (R(nx1)), u = 1x1, k = dim(phi), 
% z = 9x1 (R(k+n)), v = p = 4x1 (R(n)), r = 5x1 (dim(x)+dim(u))
% 


 
l = 0.5; m1 = 1; m2 = 0.3; g = 9.81;

% dynamics equation
f = [x3;
     x4;
    (l*m2*sin(x2)*x4^2+u1+m2*g*cos(x2)*sin(x2))/...
    (m1+m2*(sin(x2)^2));
    -(l*m2*cos(x2)*sin(x2)*x4^2+u1*cos(x2)+(m1+m2)*g*sin(x2))/...
    (l*m1+l*m2*sin(x2)^2)];
        


X_I = [x1;x2;x3;x4];                               % state vector
U_I = u1;                                          % control vector 
phi = [u1^2;x1^2;x2^2;x3^2;x4^2];                  % feature vector set


F_I = [jacobian(phi,X_I).' jacobian(f,X_I).';
     jacobian(phi,U_I).' jacobian(f,U_I).'];       % 5x9 matrix

G_I = [eye(4);zeros(1,4)];                         % 5x4 matrix


% Inverse "dynamics" matrices (will be handy to know the general shape
% of the matrices)
A_I = zeros(9);            % 9x9 matrix 
B_I = [zeros(5,4);eye(4)]; % 9x4 matrix 
 
% Inverse weight matrices
Q_I = F_I.'*F_I;   % coupling between "states"
R_I = G_I.'*G_I;   % coupling between "controls"
N_I = F_I.'*G_I;   % cross coupling



Q_fun = matlabFunction(Q_I,'vars',{x1,x2,x3,x4,u1});
N_fun = matlabFunction(N_I,'vars',{x1,x2,x3,x4,u1});



% Generating and storing time dependent Q matrices
Q_I_new = zeros(9,9,length(0:dt:2));
%Q_I_new = zeros(9,9,length(time));

for i = 1:length(0:dt:2)
%for i = 1:length(time)
    Q_I_new(:,:,i) = Q_fun(q_opt(i,1),q_opt(i,2),qdot_opt(i,1),qdot_opt(i,2),u_opt(i,1));
end

% Reverse the blocks for matrix differential Riccati equation
Q_I_bwd = Q_I_new(:,:,end:-1:1); 



% Generating and storing time dependent N matrices
N_I_new = zeros(9,4,length(0:dt:2));
%N_I_new = zeros(9,4,length(time));
for i = 1:length(0:dt:2)
%for i = 1:length(time)
    N_I_new(:,:,i) = N_fun(q_opt(i,1),q_opt(i,2),qdot_opt(i,1),qdot_opt(i,2),u_opt(i,1));
end

% Reverse the blocks for matrix differential Riccati equation
N_I_bwd = N_I_new(:,:,end:-1:1);


% Differential Ricatti Equation (with cross terms, due to residue formulation)
% P_dot = - P*A - A.'*P + P*B*inv(R)*B.'*P - Q + 0.5* P*B*inv(R)*N.' + 0.25*N*inv(R)*N.' + 0.5* N*inv(R)*B.'*P
% P(t_f) = S_f
%
% Using Euler integration to solve for matrix differential riccati equation

% initializations
P_I_new = zeros(9);



% You know P_I_new, you want to find P_I_old (backward integration)
%
% solve the differential algebraic Riccati equation


t_series = 2:-dt:0;
%t_series = flipud(time);



%[t,P] = ode45(@(t,P)matrixRiccati(t,A_I,B_I,Q_I_bwd(:,:,time_counter(t,t_series,dt)),R_I,N_I_bwd(:,:,time_counter(t,t_series,dt)),P),t_series,P_I_new);
[t,P] = ode45(@(t,P)matrixRiccati(t,A_I,B_I,Q_I_bwd,R_I,N_I_bwd,P,t_series),t_series,P_I_new);



% to pick out last five rows which is P_I(0)
% P_I_init =  P_I(end-4:end,:); 

j = length(t);
% P_mat = zeros(9,9,j);
% for i =1:j
%     P_mat(:,:,i) = reshape(P(i,:),[9,9]);
% end

% P_I_init = P_mat(:,:,end);


P_I_init = reshape(P(j,:),[9,9]);


% convex minimization
H = P_I_init;
H = (H+H.')/2;
Ainq = [0 -1 0 0 0 0 0 0 0;
        0 0 -1 0 0 0 0 0 0;
        0 0 0 -1 0 0 0 0 0;
        0 0 0 0 -1 0 0 0 0];    % size of z = [c;p]
Binq = [0;0;0;0];
Aeq  = [1 0 0 0 0 0 0 0 0];
beq = 1;

[weights,feval] = quadprog(H,[],Ainq,Binq,Aeq,beq)