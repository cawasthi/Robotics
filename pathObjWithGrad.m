function [f,g] = pathObjWithGrad(x)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Note: The difference between the full cost function and the partial cost
%       function is in the "basis" vectors (feature set). 
% 
% Full cost function    : basis :: u1,x1,x2,x3,x4
% Partial cost function : basis :: u1,x1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% 1. It seems like you cannot use a partial cost function due to the way the
%    solver works. The problem arises not in the partial cost function per 
%    se, but when you define and try to supply the corresponding partial 
%    gradient to the solver
% 2. This makes me think, will I have to reformulate the problem in case
%    I want to only do partial minimzation of the cost function.
% 3. More importantly, I need to think about how can I include new basis 
%    functions, which is probably gonna be needed when I work on the actual
%    data
%


global nGrid dt
%x1 = 4; x2 = 7; x3 = 2; x4 = 8;
x1 = 2; x2 = 3; x3 = 4; x4 = 5;


% Cost function, weighted, partial
% weights = ones(2*nGrid,1);
% weights([1,nGrid]) = 0.5;                              % quadrature weights, u1
% weights([nGrid+1,end]) = 0.5;                          % quadrature weights, x1
% weights(nGrid+1:end) = 2*weights(nGrid+1:end);         % weighted x1, cost function weights
% f= dt* x([1:nGrid,nGrid+1:4:end-3])'.^2*weights;       % objective function


% full cost function
weights = ones(5*nGrid,1);
weights([1,nGrid]) = 0.5;                                % quadrature weights, u1
weights([nGrid+1:nGrid+4,end-3:end]) = 0.5;              % quadrature weights, x1,x2,x3,x4


%weights(1:nGrid)  = 171*weights(1:nGrid);                  % TESTING WEIGHTED u1


weights(nGrid+1:4:end) = x1*weights(nGrid+1:4:end);      % weighted x1, cost function weights
weights(nGrid+2:4:end) = x2*weights(nGrid+2:4:end);      % weighted x2, cost function weights
weights(nGrid+3:4:end) = x3*weights(nGrid+3:4:end);      % weighted x3, cost function weights
weights(nGrid+4:4:end) = x4*weights(nGrid+4:4:end);      % weighted x4, cost function weights
f= dt* x'.^2*weights;                                    % objective function



% full gradient

if nargout > 1 % gradient required
    weights_u = ones(nGrid,1);
    weights_u([1,nGrid]) = 0.5;                        % quadrature weights, u1
    %weights_u(1:nGrid)   = 171*weights_u(1:nGrid);       % TESTING WEIGHTED u1
    
    
    weights_x = ones(4*nGrid,1);
    weights_x([1:4,end-3:end]) = 0.5;                  % quadrature weights, x1
    weights_x(1:4:end) = x1*weights_x(1:4:end);        % weighted x1, cost function weights
    weights_x(2:4:end) = x2*weights_x(2:4:end);        % weighted x2, cost function weights
    weights_x(3:4:end) = x3*weights_x(3:4:end);        % weighted x3, cost function weights
    weights_x(4:4:end) = x4*weights_x(4:4:end);        % weighted x4, cost function weights
    
    g = dt*2*[x(1:nGrid).*weights_u;
              x(nGrid+1:end).*weights_x];              % gradient function
               
end




% partial gradient

% if nargout > 1 % gradient required
%     weights_u = ones(nGrid,1);
%     weights_u([1,nGrid]) = 0.5;                       % quadrature weights, u1
% 
%     weights_x = ones(nGrid,1);
%     weights_x([1,nGrid]) = 0.5;                       % quadrature weights, x1
%     weights_x(1:end) = 2*weights_x(1:end);            % weighted x1, cost function weights 
% 
%     g = dt*2*[x(1:nGrid).*weights_u;
%               x(nGrid+1:4:end-3).*weights_x];         % gradient function
% end


end

