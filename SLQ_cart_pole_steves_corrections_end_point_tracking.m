function SLQ_cart_pole_steves_corrections_end_point_tracking
% Implementation of the SLQ alorithm from Fast nonlinear model predictive
% contro for unified trajectory optimization and tracking Neunert et al. ICRA 2016
%
%  Psuedocode --> page 1399 of the paper
% dimensionality check done
close all
clear  

max_slq_iter=200;
max_lin_search_iter=40;
tol=0;iter=0;

slq_tol=0.5;


% The above load x_des and u_des
 
load('energy_shaping_lqr_trajectory_cart_pole.mat');
% % dt=0.01;
% x_0=[0.0000 0.000 0.000 0.000].*1; % n_steps=1879;
% 

load('cartPole.mat');load('u_cartPole.mat'); dt=0.01; x_0=[0.0500 0.050 0.050 0.050];  n_steps=1879;

 
% states are x theta x_dot theta_dot
x_0=x_des(1,:);
x_goal=x_des(end,:);
u_des=u_des;
u_goal = 0;

Qf=1000.*eye(4); Q=1.*eye(4);  % we have 2 states so Q vectors are 2X2
R=1;                   % only 1 control so R s are scalars
% P_tf=Qf;p_tf=[1 1 1 1];
P_tf=2*Qf;
p_tf=(2*Qf*ones(4,1))';

% these are matching dimansions from the paper equation

alpha_d=1.2;             % linesearch update
u_max=8;           % control cap

% system_physical_parameters= [1;1;9.8;0.3]; % [mass; length; accln. due to gravity; damping]

iter=0;
del_l=999;
del_u=999;
del_J=999;
lt=9999.*ones(n_steps,1); %  || del_l >=  slq_tol ||
J_now=999;
figure(1)
xlabel('Pole angle');
ylabel('Pole ang. velocity');
hold all

% making the desired trajectories same as being at goal at all times
x_nom=x_des;
u_nom=u_des ;
  plot(x_nom(:,2),x_nom(:,4),'r.-');
x_des=repmat(x_goal,[max(size(x_des)),1]);
u_des=zeros(max(size(u_des)),1);
 
while  iter< max_slq_iter  &&  del_l >=  slq_tol
    %     disp('new iter');
    l_old=lt;
    u_old=u_nom;
    % this shows the status at each iteration
    x_nom=simulate_dynamics(u_nom, x_0,dt); %we do not need a x_nom tape as we generate it here by rollout
    plot(x_nom(:,2),x_nom(:,4),'b.-');
    pause
    %         hold all
    % to calculate updates for r, q
    %     initialize
    sq=1.*ones(n_steps,1);
    r=ones(n_steps,1);
    for i=1:n_steps
        %         sq(i,:)=-(2*Q*(x_des(i,:)-x_nom(i,:))')';
        %         r(i)     =-dot(2*diag(R),(u_des(i)-u_nom(i)));
        q(i,:)=(2*Q*Wrap_difference(x_nom(i,:),x_des(i,:),2)')';
        r(i)=dot(2*diag(R),(u_nom(i)-u_des(i)));
    end
    % or change q and r to ones by commenting the for loop!
    %     q=sq ;%ones(n_steps,4); %r=ones(n_steps,1);
    
    [Kt,lt]= riccati_like(P_tf, p_tf,Q,q,R,r,x_des,x_nom,u_nom,dt);
    
    J_now=total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des);
    
    J_new=inf;
    u_iter=1;
    alpha=1;
    J_buff=[0];
    u_buff=zeros(1,n_steps);
    J_buff(u_iter)=J_new;
    J_u_buff(u_iter,:)=u_nom(:)';
    %      && J_new  > J_now
    x_0_ct= x_0;
    while  u_iter <= max_lin_search_iter && J_new  >= J_now
        %     update control step by step
%         for i=1:n_steps
%             %             u_nom(i)=u_old(i)+alpha*lt(i) +
%             %             Kt(i,:)*(x_des(i,:)-x_nom(i,:))'; Wrap_difference(x_0_ct , x_nom(i,:),2)'
%             % %             u_nom(i)=u_old(i) + alpha*lt(i) +  Kt(i,:)*(x_0_ct - x_nom(i,:))';  % to do wrapto pi for angle
%             u_nom(i)=u_old(i) + alpha*lt(i) +  Kt(i,:)*Wrap_difference(x_0_ct, x_des(i,:),2)';
%             if u_nom(i)>=u_max
%                 u_nom(i)=u_max;
%             elseif u_nom(i)<=-u_max
%                 u_nom(i)=-u_max;
%             end
%             x_0_ct = f(x_0_ct,u_nom(i))'.*dt + x_0_ct;
%         end

%  update control all at once
        for i=1:n_steps
            u_old(i) + alpha*lt(i) +  Kt(i,:)*Wrap_difference(x_nom(i,:), x_des(i,:),2)';
            if u_nom(i)>=u_max
                u_nom(i)=u_max;
            elseif u_nom(i)<=-u_max
                u_nom(i)=-u_max;
            end
            
        end
        
        %-- control update done--
        
        x_nom=simulate_dynamics(u_nom, x_0,dt);
        J_new=total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des);
        J_buff(u_iter+1)=J_new;
        u_buff(u_iter+1,:)=u_nom(:)';
        alpha=alpha/alpha_d;
        u_iter=u_iter+1;
    end
    
    %     size(u_buff)
    % we choose the u corresponding to the minimum J_now
    u_iter
    [~,min_id]=min(J_buff)
    u_nom=u_buff(min_id,:);
    J_new=J_buff(min_id);
    
    %     Iteration of the outer while loop
    iter=iter+1
    %     del_J=J_new-J_now;
    del_l=max(abs(l_old-lt))  % our termination condition
    %     del_u=norm(u_old-u_nom,2)
    
end


%  plot to show how we did
hold on
x_nom=simulate_dynamics(u_nom, x_0,dt); % rollout with the final solution
plot(x_nom(:,2),x_nom(:,4),'b.-');
plot(x_des(:,2),x_des(:,4),'r.-');
figure(2)
plot(u_nom)
title('control effort');
figure(3)
title('cart phase plot')
plot(x_nom(:,1),x_nom(:,3),'b.-');
xlabel('cart position');
ylabel('cart velocity');

figure(4)
x_nom=[x_nom;x_nom(end,:)];
plot_cart_pole (x_nom',dt.*[1:1:n_steps+1]);
end


function [Kt_out,lt_out]= riccati_like(P_tf, p_tf,Q,q,R,r,x_des,x_nom,u_nom,dt)
%  first populate the variable vectors
% Q  and q are the full state cost
% R and r are the control costs
% P is the quadratic difference cost
% p is the linear difference cost
[n_steps,~]=size(x_nom); %2 x vector
P=zeros(n_steps,16); % we will reshape this to [2,2] before use
P(end,:)=P_tf(:)'; % from the input
p=zeros(n_steps,4); %will reshape this into [2,1]
p(end,:)=p_tf(:)';% from the input
m=1;l=1;g=1;b=1; % system properties

% initialize outputs
Kt_out=zeros(n_steps,4);
lt_out=zeros(n_steps,1);

t_vec=[2:1:n_steps];

for i = fliplr(t_vec)
    u=u_nom(i);
    [A_continuous,B_continuous] = AB_Partial(x_nom(i,:),u);
    
    %     M=[A_continuous B_continuous; zeros(1,5)].*dt;
    %     Mt=expm(M);
    %     At=Mt(1:4,1:4);
    %     Bt=Mt(1:4,5);
    [A_discrete,B_discrete] = linearize(A_continuous,B_continuous,dt) ;
    At = A_discrete ;
    Bt = B_discrete ;
    
    
    
    %     [Atstr,Btstr] = Discrete_A_B(x_nom(i,:),u_nom(i),dt);
    %     Atstr-At
    %     Btstr-Bt
    
    P_t_plus_1 = reshape(P(i,:),[4,4]); %% 2X2 square matrix
    p_t_plus_1 = p(i,:)';  %% 2X1 vector
    
    Ht=R+Bt'*P_t_plus_1 *Bt; %% scalar-- 1+ 1X2*2X2*2X1
    Gt=Bt'*P_t_plus_1*At ; %% 1X2 vector --1X2*2X2*2X2
    gt=r(i)+Bt'*p_t_plus_1;   %% scalar -- 1+1X2*2X1
    
    Kt=-inv(Ht)*Gt ; %% 1X2 vector--  1*1X2
    Kt_out(i,:)=Kt(:); %% for storage!-- we need to reshape them
    %     pause
    
    %     lt_out(i,:)= - inv(Ht)*gt - inv(Ht)*Gt*(x_nom(i,:)-x_des(i,:))';  %% scalar 1*1 --> 1*1 -1* 1X2*2X1
    lt_out(i,:)= - inv(Ht)*gt - inv(Ht)*Gt*Wrap_difference(x_nom(i,:),x_des(i,:),2)';
    lt=lt_out(i,:);           %% for storage!
    
    P_sq=Q + At'*P_t_plus_1*At + Kt'*Ht'*Kt + Kt'*Gt + Gt'*Kt;
    % % 2X2 + 2X2*2X2*2X2 + 2X1*1*1X2  +  2X1*1X2  +  2X1*1X2  = 2X2
    P(i-1,:)=P_sq(:)'; %% for storage!-- we need to reshape them -- done in line 68
    
    p_sq=q(i,:)' + At'*p_t_plus_1 + Kt'*Ht*lt + Kt'*gt + Gt'*lt; % Namam's correction
    % % 2X1 = 2X1 + 2X2*2X1 + 2X1*1*1  + 2X1*1 + 2X1*1
    
    %     p_sq=q+At'*p_t_plus_1+Kt'*Ht*lt+lt'*gt+Gt'*lt;  % original paper equation
    
    p(i-1,:)=p_sq';
end
end

%%
function [Sim_x]=simulate_dynamics(u_in, x_0,dt)
% simulate the dunamics for the damped simple pendulum
n_steps=length(u_in);
Sim_x=zeros(n_steps,4); % as we have 2 states here
Sim_x(1,:)=x_0;
for i=1:n_steps-1
    Sim_x(i+1,:)=f(Sim_x(i,:),u_in(i))'.*dt +Sim_x(i,:);
    Sim_x(i+1,2)=sign(Sim_x(i+1,2))*wrapTo2Pi(abs(Sim_x(i+1,2)));
end

end

%%
function state_dot = f(state,u)
% m=system_physical_dynamics(1);
% l=system_physical_dynamics(2);
% g=system_physical_dynamics(3);
% b=system_physical_dynamics(4);
% x=state(1);theta=state(2); x_dot=state(3); theta_dot=state(4);
% m=1;l=1;g=1;b=1; %% these are the parameters
% % At= [0 1; -g/l*cos(x(1)) -b];
% % Bt=[0 1];

ss = sin(state(2));
cc = cos(state(2));
xddot = (u+ss*state(4)^2+ss*cc)/(2-cc*cc);
tddot = (-u*cc - state(4)^2*cc*ss - 2*ss)/(2-cc*cc);
state_dot = [state(3);state(4); xddot; tddot];


end

%%
function J= total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des)

del_x_f=(x_des(end,:)-x_nom(end,:));

h = del_x_f*p_tf' + 0.5.* del_x_f*P_tf*del_x_f' ;%terminal stuff
cost=0;

for i=1:numel(u_nom) % cost to go
%         del_x=(x_nom(i,:)-x_des(i,:))'
    del_x=Wrap_difference(x_nom(i,:),x_des(i,:),2)';
%     pause
    del_u=(u_nom(i)-u_des(i));
    cost=cost + sq(i,:) + del_x'*q(i,:)' + del_u'*r(i) + 0.5.*del_x'*Q*del_x + del_u'*R*del_u;
end
% cost;
J=h+cost; % total cost
% pause
end

%%
% auxiliary functions
function [At,Bt] = Discrete_A_B(x_nom,u,dt)
% This is not used. This is just a symbolic version of the 2 term Taylor
% series expansion. We used this to see how the discrete system works!
x=x_nom(1);theta=x_nom(2); x_dot=x_nom(3); theta_dot=x_nom(4);

At= [1 0 dt 0;
    0 1 0 dt;
    0 dt * ((cos(theta) * theta_dot ^ 2 - sin(theta) ^ 2 + cos(theta) ^ 2) / (0.1e1 + sin(theta) ^ 2) - 0.2e1 * (sin(theta) * theta_dot ^ 2 + u + cos(theta) * sin(theta)) / (0.1e1 + sin(theta) ^ 2) ^ 2 * cos(theta) * sin(theta)) 1 0.2e1 * dt * sin(theta) * theta_dot / (0.1e1 + sin(theta) ^ 2);
    0 dt * ((sin(theta) ^ 2 * theta_dot ^ 2 - cos(theta) ^ 2 * theta_dot ^ 2 + u * sin(theta) - 0.2e1 * cos(theta)) / (0.1e1 + sin(theta) ^ 2) - 0.2e1 * (-cos(theta) * sin(theta) * theta_dot ^ 2 - u * cos(theta) - 0.2e1 * sin(theta)) / (0.1e1 + sin(theta) ^ 2) ^ 2 * cos(theta) * sin(theta)) 0 0.1e1 - 0.2e1 * dt * cos(theta) * sin(theta) * theta_dot / (0.1e1 + sin(theta) ^ 2);];

Bt=[0;0;0;0];


Bt(1)= dt ^ 3 / (0.1e1 + sin(theta) ^ 2) / 0.2e1;

Bt(2)= -dt ^ 3 * cos(theta) / (0.1e1 + sin(theta) ^ 2) / 0.2e1;

Bt(3)= dt ^ 2 / (0.1e1 + sin(theta) ^ 2) / 0.2e1 - dt ^ 3 * sin(theta) * theta_dot / (0.1e1 + sin(theta) ^ 2) ^ 2 * cos(theta);

Bt(4)=  -(0.1e1 - 0.2e1 * dt * cos(theta) * sin(theta) * theta_dot / (0.1e1 + sin(theta) ^ 2)) * dt ^ 2 * cos(theta) / (0.1e1 + sin(theta) ^ 2) / 0.2e1;

end


%%
function [A_continuous,B_continuous] = AB_Partial(x,u)
% Thanks Steven Crews!
A_continuous = zeros(4,4);
B_continuous = zeros(4,1);
eps = 1e-4;

g=f(x,u) ;

for i = 1:4
    x_ = x ;
    x_(i) = x_(i)+eps ;
    g_eps=f(x_,u) ;
    A_continuous(:,i)=(g_eps-g)/eps ;
end

for i = 1:1
    u_ = u ;
    u_(i) = u_(i)+eps ;
    g_eps=f(x,u_) ;
    B_continuous(:,i)=(g_eps-g)/eps ;
end
end

%%
% discretization of A and B matrices
function [A_discrete,B_discrete] = linearize(A_continous,B_continous,dt)
%  Steve Crews, code adopted from Saumya Saxena
%  Biorobotics Lab
%  Carnegie Mellon University

Nx = 4 ;% number of states
Nu = 1 ;% number of controls

M = [A_continous B_continous; zeros(1,Nx+Nu)].*dt ;
MM = expm(M) ;
A_discrete = MM(1:Nx,1:Nx);
B_discrete = MM(1:Nx,Nx+1:end);
end

function wrapped_state = Wrap_difference(x_nom,x_des,n_wrap)
wrapped_state=zeros(1,numel(x_des));
for i=1:numel(x_des)
    if i~=n_wrap
        wrapped_state(i)=x_nom(i)-x_des(i);
    else
        wrapped_state(i)=wrapTo2Pi(x_nom(i)-x_des(i));
    end
end

end
