function SLQ_MPC
close all
clear all
%  Given system dynamics

% set parameters
% 1. load target trajectory generated from dircol
load('data_for_Q4.mat');
dt=0.025; 
x_des=[theta' theta_dot'];
u_des=2*optimal;
plot(x_des(:,1),x_des(:,2),'r.-');
hold all

n_steps=numel(optimal); %% total number of timesteps through which the system is tracked

% 2. Set MPC params
n_sim=3; % roll out current strategy through so many steps through time
x_start= x_des(1,:);
x_start=[0.5 0.0];%x_start=x0; % this is where we start from. Should be in the neighborhood of 
% start point of the known optimal path. So that the lqr can attract it to
% the start of the optimal trajectory
x_goal=[pi,0];
iter=0; % required to keep track of time horizon
Q_lqr= 1500.*eye(2);
Q_slq=1500.*eye(2);
R_lqr=1; 
R_slq=1;
% these 2 parameters will be used for the LRQ ricatti stuff

x_goal=[pi,0]; % always fixed for the goal state
% begin MPC loop
%  [u_slq]= SLQ_solve(x_des,u_des,x_start,Q_slq,R_slq,100.*eye(2),[1:1:200],[2 2],dt); % to see if slq works
%  hold all
%  pause
% intermediate_counter
while (1) % we will run this until finished
    [At,Bt] = lin_dynamics(x_start);% linear dynamics at start point of the known optimal path
    [K,~,~] = lqr(At,Bt,Q_lqr,R_lqr);
    t_left=n_steps-(iter+n_sim);
    iter=iter+1;
    if t_left>=n_sim
        t_vec=[iter:1:iter+n_sim];
    else
%         intermediate_counter=intermediate_counter+1;
        t_vec=[iter:1:n_steps]
        pause
    end
t_left
%   we need to initialize the SLQ with x_start and the inf horizon gains from LQR    
    [At,Bt] = lin_dynamics(x_des(t_vec(end),:));
    [~,Q_f_SLQ,~] = lqr(At,Bt,Q_lqr,R_lqr); % the final cost from ricatti
    
    [u_slq]= SLQ_solve(x_des,u_des,x_start,Q_slq,R_slq,Q_f_SLQ,t_vec,K,dt);
%     figure(2)
%     plot(u_slq)
%     send control to forward rollout
    u_size=length(u_slq);
    if u_size>= n_sim
%         pause
        u_rollout=u_slq(1:n_sim);
    else
        disp('error!!')
        u_rollout=u_slq;
    end
    
%    u_rollout
    [Sim_x]=simulate_dynamics(u_rollout(1),x_start,dt);
    x_start=Sim_x(end,:);
    size(Sim_x)
    plot(Sim_x(:,1),Sim_x(:,2),'g.-');
    hold all
    
pause
    if numel(t_vec) <=1 % no more time steps left
        disp('Finished');
        return
    end
end

end

function [At,Bt] = lin_dynamics(x)
% m=system_physical_dynamics(1);
% l=system_physical_dynamics(2);
% g=system_physical_dynamics(3);
% b=system_physical_dynamics(4);

m=1;l=1;g=1;b=1;
At= [0 1; -g/l*cos(x(1)) -b];
Bt=[0; 1];

end

function [Sim_x]=simulate_dynamics(u_in, x_0,dt)
% simulate the dunamics for the damped simple pendulum
n_steps=length(u_in);

Sim_x=zeros(n_steps,2); % as we have 2 states here
Sim_x(1,:)=x_0;
for i=1:n_steps
    Sim_x(i+1,:)=f(Sim_x(i,:),u_in(i))'.*dt +Sim_x(i,:);
%             Sim_x(i+1,1)=wrapTo2Pi(Sim_x(i+1,1));
end

end
function x_dot = f(x,u)
% m=system_physical_dynamics(1);
% l=system_physical_dynamics(2);
% g=system_physical_dynamics(3);
% b=system_physical_dynamics(4);

m=1;l=1;g=1;b=1;
% % At= [0 1; -g/l*cos(x(1)) -b];
% % Bt=[0 1];
x_dot=[x(2);-g*sin(x(1))/l - b*x(2)/(m*l*l) + u/(m*l^2)];
end