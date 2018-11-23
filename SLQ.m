function SLQ
% Implementation of the SLQ alorithm from Fast nonlinear model predictive
% contro for unified trajectory optimization and tracking Neunert et al. ICRA 2016
%
%  Psuedocode --> page 1399 of the paper
% dimensionality check done
close all
clear all
load('data_for_Q4.mat');
max_slq_iter=10;
max_lin_search_iter=10;
tol=0;iter=0;
slq_tol=0.001;
cost_tol=0.001; tol_cost=999;
n_states=200;
x_des=[theta_sol' theta_dot_sol'];
u_des=u0;
x_0=[0 0]; % also equal to the first row of the x_nom
x_goal=[pi,0];
Qf=eye(2); Q=eye(2);q=[1;0];  % we have 2 states so Q vectors are 2X2
R=1;r=1;                   % only 1 control so R s are scalars
P_tf=eye(2);p_tf=[1 1];   % these are matiching dimansions from the paper equation
alpha_d=0.1;             % linesearch update
% system_physical_parameters= [1;1;9.8;0.3]; % [mass; length; accln. due to gravity; damping]
iter=0;
while  iter< max_slq_iter
    x_nom=simulate_dynamics(u_nom, x_0);
    figure(1)
    hold all
    % size(x_nom)
    plot(x_nom(:,1),x_nom(:,2),'r.-');
    [Kt,lt]= riccati_like(P_tf,p_tf,x_nom,Q,q,R,r);
    alpha=1;
    J_now=total_cost(x_nom, u_nom, Qf,Q,R,x_goal)
    pause
    J_new=1e10000;
    u_iter=0;
    while J_new >= J_now
        %     update control
        u_new=zeros(1,n_states);
        for i=1:n_states
            u_new(i)=u_nom(i)+alpha*lt(i)+Kt(i,:)*(x_des(i,:)-x_nom(i,:))';
        end
        x_nom=simulate_dynamics(u_new, x_0);
        J_new=total_cost(x_nom, u_nom, Qf,Q,R,x_goal);
        alpha=alpha/alpha_d;
        u_iter=u_iter+1;
    end
    iter=iter+1
    pause
end

end


function [Kt_out,lt_out]= riccati_like(P_tf, p_tf,x_tape,Q,q,R,r)
%  first populate the variable vectors
% Q  and q are the full state cost
% R and r are the control costs
% P is the quadratic difference cost
% p is the linear difference cost
[n_steps,~]=size(x_tape); %2 x vector
P=zeros(n_steps,4); % we will reshape this to [2,2] beefore use
P(end,:)=[P_tf(1,1) P_tf(1,2) P_tf(2,1) P_tf(2,2)]; % from the input
p=zeros(n_steps,2); %will reshape this into [2,1]
p(end,:)=[p_tf(1) P_tf(2)];% from the input
m=1;l=1;g=9.8;b=0.3; % system properties

% initialize outputs
Kt_out=zeros(n_steps,2);lt_out=zeros(n_steps,1);

t_vec=[2:1:n_steps];
for i = fliplr(t_vec)
    Bt=[0;1]; %% 2X1 vector
    At= [0 1; -g/l*cos(x_tape(i,1)) -b]; %% 2X2 square matrix
    P_t_plus_1 = reshape(P(i,:),[2,2]); %% 2X2 square matrix
    p_t_plus_1 = p(i,:)';  %% 2X1 vector
    
    Ht=R+Bt'*P_t_plus_1 *Bt; %% scalar-- 1+ 1X2*2X2*2X1
    Gt=Bt'*P_t_plus_1*At ; %% 1X2 vector --1X2*2X2*2X2
    gt=r+Bt'*p_t_plus_1;   %% scalar -- 1+1X2*2X1
    
    
    Kt=-inv(Ht)*Gt';  %% 2X1 vector--  2X2*2X1
    Kt_out(i,:)=Kt(:)'; %% for storage!-- we need to reshape them
    
    lt_out(i,:)=-inv(Ht)*gt;  %% scalar 1*1
    lt=lt_out(i,:);           %% for storage!
    
    P_sq=Q + At'*P_t_plus_1*At + Kt*Ht'*Kt' + Kt*Gt + Gt'*Kt';
    % % 2X2 + 2X2*2X2*2X2 + 2X1*1*1X2  +  2X1*1X2  +  2X1*1X2  = 2X2
    P(i-1,:)=P_sq(:)'; %% for storage!-- we need to reshape them -- done in line 68
    
    p_sq=q + At'*p_t_plus_1 + Kt*Ht*lt + Kt*gt + Gt'*lt; % Namam's correction
    % % 2X1 = 2X1 + 2X2*2X1 + 2X1*1*1  + 2X1*1 + 2X1*1
    
    %     p_sq=q+At'*p_t_plus_1+Kt'*Ht*lt+lt'*gt+Gt'*lt;  % original paper equation
    
    p(i-1,:)=p_sq';
end
end
function [Sim_x]=simulate_dynamics(u_in, x_0)
% simulate the dunamics for the damped simple pendulum
n_steps=length(u_in);
n_states=length(x_0);
Sim_x=zeros(n_states,2); % as we have 2 states here
Sim_x(1,:)=x_0;
for i=2:n_steps
    State_prev=Sim_x(i-1,:);
    Sim_x(i,:)=f(State_prev,u_in(i));
end
end
function x_dot = f(x,u)
% m=system_physical_dynamics(1);
% l=system_physical_dynamics(2);
% g=system_physical_dynamics(3);
% b=system_physical_dynamics(4);
m=1;l=1;g=9.8;b=0.3;
At= [0 1; -g/l*cos(x(1)) -b];
Bt=[0 1];
x_dot=At*x'+Bt'*u;
end

function J= total_cost(x,u,Qf,Q,R,P_tf,p_tf,x_des,u_des)
del_x_f=x(end,:)-x_des(end,:)
del
h=(x(end,:)-x_goal)*Qf*(x(end,:)-x_goal)'
cost=0;
for i=1:numel(u)
    cost=cost+x(i,:)*Q*x(i,:)'+u(i)*R*u(i)';
end
J=h+cost;
end
