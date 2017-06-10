function [c ,ceq] = dynConsG(x) 

global nGrid dt
c = [];


% % Getting states/control out of the vector x
% Getting states/control out of the vector x
u  = x(1:nGrid);            % 30x1 (10x3)
q1 = x(nGrid+1:4:end-3);    % 10x1
q2 = x(nGrid+2:4:end-2);    % " "     
q3 = x(nGrid+3:4:end-1);    % " "          
q4 = x(nGrid+4:4:end);      % " "     


q = [q1,q2]';      % 3x10
q_dot = [q3,q4]';  % 3x10


% NOTES: lot of confusiing terminology floating around with
% q_dot, f_dot,x_dot
% q_dot = joint velocity
% f_dot = "state" derivative
% x_dot = overall "state" derivative (which gets passed to fmincon), I don't think we need x_dot
% f_i should be used instead of x_i below    
    
    
% Constrain initial position and velocity to be zero
ceq = [q(:,1);q_dot(:,1)];


for i = 1 : nGrid - 1
    
    % The state at the beginning of the time interval
    f_i = [q(:,i); q_dot(:,i)];   % 6x1
    
    
    % What the state should be at the start of the next time interval
    f_n = [q(:,i+1); q_dot(:,i+1)];
    
    
    % The time derivative of the state at the beginning of the time
    % interval
    fdot_i =  robDyn(q(:,i),q_dot(:,i),u(i)); 
         
        
    % The time derivative of the state at the end of the time interval
    fdot_n = robDyn(q(:,i+1),q_dot(:,i+1),u(i+1));  
    
        
    % The end state of the time interval calculated using quadrature
    f_end = f_i + dt * (fdot_i + fdot_n) / 2;
    
    
    % Constrain the end state of the current time interval to be
    % equal to the starting state of the next time interval
    
    ceq = [ceq ; f_n - f_end];          %#ok<AGROW>
    
  
end

% Constrain end position to given and end velocity to zero
ceq = [ceq ; q(:,end) - [1;pi]; q_dot(:,end)];

  
end  




function [f_dot] = robDyn(q,q_dot,u)


% system dynamics
l = 0.5; m1 = 1; m2 = 0.3; g = 9.81;

% dynamics equation
f_Dot = [q_dot(1);
        q_dot(2);
        (l*m2*sin(q(2))*q_dot(2)^2+u+m2*g*cos(q(2))*sin(q(2)))/...
        (m1+m2*(sin(q(2))^2));
        -(l*m2*cos(q(2))*sin(q(2))*q_dot(2)^2+u*cos(q(2))+(m1+m2)*g*sin(q(2)))/...
        (l*m1+l*m2*sin(q(2))^2)];
        

    
q_dot = f_Dot(1:2,1);
q_ddot = f_Dot(3:4,1);


f_dot = [q_dot ;q_ddot];

end

