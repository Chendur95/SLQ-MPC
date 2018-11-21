function SLQ
% Implementation of the SLQ alorithm from Fast nonlinear model predictive
% contro for unified trajectory optimization and tracking Neunert et al. ICRA 2016 
% 
%  Psuedocode --> page 1399 of the paper

max_slq_iter=10;
max_lin_search_iter=10;
slq_tol=0.001;
lin_search_tol=0.001;
% system_physical_parameters= [1;1;9.8;0.3]; % [mass; length; accln. due to gravity; damping]

[simulated_states]=simulate_dynamics(u_in,system_physical_params); % linearised system dymanics inside

end
function [P,K,l]= riccati_like(P_tf, p_tf,x_tape,Q,q,R,r)
%  first populate the variable vectors
% Q  and q are the full state cost
% R and r are the control costs
% P is the quadratic difference cost
% p is the linear difference cost
[~,n_steps]=size(x_tape); %2 x vector 
P=zeros(n_steps,4); % we will reshape this to [2,2] beefore use
P(end,:)=[P_tf(1,1) P_tf(1,2) P_tf(2,1) P_tf(2,2)]; % from the input
p=zeros(n_steps,2); %will reshape this into [2,1]
p(end,:)=[p_tf(1) P_tf(2)];% from the input
m=1;l=1;g=9.8;b=0.3; % system properties

for i=n_steps:1
    Bt=[0;1];
    At= [0 1; -g/l*cos(x_tape(i,1)) -b];
    P_t_plus_1 = reshape(P(i,:),[2,2]);
    p_t_plus_1 = p(i,:)';
    
    Ht=R+Bt'*P_t_plus_1 *Bt;
    Gt=Bt'*P_t_plus_1*At;
    gt=r+Bt'*p_t_plus_1;
    
    Kt
    
    P_sq=Q+At'*P_t_plus_1*At+Kt'*Ht'*Kt+Kt'*Gt+Gt'*Kt; 
    P(i,:)=P_sq(:)';
    
    p_sq=q+At'*p_t_plus_1+Kt'*Ht*lt+lt'gt+Gt'lt;
end
end
function [simulated_states]=simulate_dynamics(u_in, x_0)
% simulate the dunamics for the damped simple pendulum
n_steps=length(u_in);
nstates=length(x_0);
Sim_x=zeros(n+1,n_states); % as we have 2 states here
Sim_x(1,:)=x_0;
for i=2:n
%    del_x= Sim_x(i)-Sim_x(i-1);
%    del_u=u_in(i)-u_in(i-1);
   Sim_x(i,:)=Sim_x(i-1,:)+f(Sim_x(i-1,:),del_x,del_u)';
end
end
function x_dot = f(x_t,x,u)
% m=system_physical_dynamics(1);
% l=system_physical_dynamics(2);
% g=system_physical_dynamics(3);
% b=system_physical_dynamics(4);
m=1;l=1;g=9.8;b=0.3;
At= [0 1; -g/l*cos(x_t(1)) -b];
Bt=[0 1];
x_dot=At*x+Bt*u;
end
