function SLQ_cart_pole
% Implementation of the SLQ alorithm from Fast nonlinear model predictive
% contro for unified trajectory optimization and tracking Neunert et al. ICRA 2016
%
%  Psuedocode --> page 1399 of the paper
% dimensionality check done
close all
clear all

max_slq_iter=50;
max_lin_search_iter=20;
tol=0;iter=0;

slq_tol=1;
cost_tol=8;

%  input is just the control tape
% load('dir_col_trajectory_cart_pole.mat'); u_des=u.*0;
% load('energy_shaping_lqr_trajectory_cart_pole.mat');n_steps=800;
load('cartPole.mat');load('u_cartPole.mat'); dt=0.01; x_0=[0.0500 0.050 0.050 0.050];  n_steps=1879;

% dt loaded from file
% u_loaded from file
% x_des loaded from file
% n_steps loaded from file
% x_des=zeros(800,4);
% u_des=zeros(800,1);
% % n_steps=800;
% n_steps=1879;
% % dt=0.05;
%  u_nom=u_des;
%  dt=0.01;
% x_des loaded from file
% hold all
% u_nom=0.001.*ones(n_steps,1);
% u_des=zeros(500,1);
% % also equal to the first row of the x_nom
% states are x theta x_dot theta_dot
% for ii=2:numel(u_des)
%     [At,Bt] = Discrete_A_B(x_0(ii-1,:),u_des(ii-1),dt);
%     x_0(ii,:)=At*x_0(ii-1,:)'+Bt*u_des(ii-1);
% end
% plot(x_0(:,1),'g.-');pause
% figure (1)
% plot(x_des(:,2),x_des(:,4),'g.-');
% hold all

% x_str=simulate_dynamics(u_des, x_0,dt);
% % size(u_des)
% % size(x_str)
% plot(x_str(:,2),x_str(:,4),'r.-');

% pause
% x_0=[0.000 0.000 0.000 0.000];
x_goal=[0,pi,0,0];
Qf=1.*eye(4); Q=1.*eye(4);  % we have 2 states so Q vectors are 2X2
R=1;                   % only 1 control so R s are scalars
P_tf=Qf;p_tf=[1 1 1 1];
% these are matching dimansions from the paper equation
alpha_d=0.5;             % linesearch update
u_max=10;

% system_physical_parameters= [1;1;9.8;0.3]; % [mass; length; accln. due to gravity; damping]
u_nom=u_des*1 ;
iter=0;
del_l=999;
del_u=999;
del_J=999;
lt=0.*ones(n_steps,1); %  || del_l >=  slq_tol ||
J_now=999;
while  iter< max_slq_iter  &&  del_l >=  slq_tol && J_now>= slq_tol
    disp('new iter');
    l_old=lt;
    u_old=u_nom;
    x_nom=simulate_dynamics(u_nom, x_0,dt); %we do not need a x_nom tape as we generate it here by rollout
    plot(x_nom(:,2),x_nom(:,4),'b.-');
    pause
    %         hold all
    % to calculate updates for r, q
    %     initialize
    sq=1.*ones(n_steps,1); r=ones(n_steps,1);
            for i=1:n_steps
                sq(i,:)=2*diag(Q)'*(x_des(i,:)-x_nom(i,:))';
                r(i)=dot(2*diag(R),(u_des(i)-u_nom(i)));
            end
    % or change q and r to ones!
    q=ones(n_steps,4); %r=ones(n_steps,1);
    
    
    [Kt,lt]= riccati_like(P_tf, p_tf,Q,q,R,r,x_des,x_nom,u_nom,dt);
    
    J_now=total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des)
    
    J_new=inf;
    u_iter=1;  %  J_new -J_now >= cost_tol | && J_new  >= J_now ||
    alpha=1;
    J_buff=[0];
    u_buff=zeros(1,n_steps);
    J_buff(u_iter)=J_new;
    u_buff(u_iter,:)=u_nom(:)';
    
    while  u_iter <= max_lin_search_iter && J_new  > J_now
        %     update control
        for i=1:n_steps
            u_nom(i)=u_old(i)+alpha*lt(i) + Kt(i,:)*(x_des(i,:)-x_nom(i,:))';
            if u_nom(i)>=u_max
                u_nom(i)=u_max;
            elseif u_nom(i)<=-u_max
                u_nom(i)=-u_max;
            end
            
        end
        
        x_nom=simulate_dynamics(u_nom, x_0,dt);
        J_new=total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des);
        J_buff(u_iter+1)=J_new;
        u_buff(u_iter+1,:)=u_nom(:)';
        alpha=alpha/alpha_d;
        u_iter=u_iter+1;
    end
    
%     size(u_buff)
    [~,min_id]=min(J_buff)
    u_nom=u_buff(min_id,:);
    J_new=J_buff(min_id);
    
    iter=iter+1
%     del_J=J_new-J_now;
    del_l=max(abs(l_old-lt))
%     del_u=norm(u_old-u_nom,2)
    
end
%  plot to show how we did
hold on
x_nom=simulate_dynamics(u_nom, x_0,dt);
plot(x_nom(:,2),x_nom(:,4),'b.-');
plot(x_des(:,2),x_des(:,4),'r.-');
figure(2)
plot(u_nom)
figure(3)
plot(x_nom(:,1),x_nom(:,3),'b.-');
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
    
%     x=x_nom(i,1);theta=x_nom(i,2); x_dot=x_nom(i,3); theta_dot=x_nom(i,4);
%     A=[0 0 1 0; 0 0 0 1; 0 (cos(theta) * theta_dot ^ 2 - sin(theta) ^ 2 + cos(theta) ^ 2) / (0.1e1 + sin(theta) ^ 2) - (0.2e1 * sin(theta) * theta_dot ^ 2 + (2 * u) + 0.2e1 * cos(theta) * sin(theta)) * cos(theta) * sin(theta) / (0.1e1 + sin(theta) ^ 2) ^ 2 0 0.2e1 * sin(theta) * theta_dot / (0.1e1 + sin(theta) ^ 2); 0 -(-sin(theta) ^ 2 * theta_dot ^ 2 + cos(theta) ^ 2 * theta_dot ^ 2 - u * sin(theta) + 0.2e1 * cos(theta)) / (0.1e1 + sin(theta) ^ 2) + (0.2e1 * cos(theta) * sin(theta) * theta_dot ^ 2 + 0.2e1 * u * cos(theta) + 0.4e1 * sin(theta)) * cos(theta) * sin(theta) / (0.1e1 + sin(theta) ^ 2) ^ 2 0 -0.2e1 * cos(theta) * sin(theta) * theta_dot / (0.1e1 + sin(theta) ^ 2);];
%     B=[0; 0; 0.1e1 / (0.1e1 + sin(theta) ^ 2); -cos(theta) / (0.1e1 + sin(theta) ^ 2)];
        
[A,B] = AB_Partial(x_nom(i,:),u);
    
    M=[A B; zeros(1,5)].*dt; Mt=expm(M); At=Mt(1:4,1:4); Bt=Mt(1:4,5);
    
    [Atstr,Btstr] = Discrete_A_B(x_nom(i,:),u_nom(i),dt);
    
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
    
    lt_out(i,:)= - inv(Ht)*gt - inv(Ht)*Gt*(x_des(i,:)-x_nom(i,:))';  %% scalar 1*1 --> 1*1 -1* 1X2*2X1
    lt=lt_out(i,:);           %% for storage!
    
    P_sq=Q + At'*P_t_plus_1*At + Kt'*Ht'*Kt + Kt'*Gt + Gt'*Kt;
    % % 2X2 + 2X2*2X2*2X2 + 2X1*1*1X2  +  2X1*1X2  +  2X1*1X2  = 2X2
    P(i-1,:)=P_sq(:)'; %% for storage!-- we need to reshape them -- done in line 68
    
    p_sq=q(i) + At'*p_t_plus_1 + Kt'*Ht*lt + Kt'*gt + Gt'*lt; % Namam's correction
    % % 2X1 = 2X1 + 2X2*2X1 + 2X1*1*1  + 2X1*1 + 2X1*1
    
    %     p_sq=q+At'*p_t_plus_1+Kt'*Ht*lt+lt'*gt+Gt'*lt;  % original paper equation
    
    p(i-1,:)=p_sq';
end
end
function [Sim_x]=simulate_dynamics(u_in, x_0,dt)
% simulate the dunamics for the damped simple pendulum
n_steps=length(u_in);

Sim_x=zeros(n_steps,4); % as we have 2 states here
Sim_x(1,:)=x_0;
for i=1:n_steps-1
    Sim_x(i+1,:)=f(Sim_x(i,:),u_in(i))'.*dt +Sim_x(i,:);
       Sim_x(i+1,1)=wrapTo2Pi(Sim_x(i+1,1));
end

end
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

function J= total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des)

del_x_f=(x_des(end,:)-x_nom(end,:));

h = del_x_f*p_tf' + 0.5.* del_x_f*P_tf*del_x_f' ;%terminal stuff
cost=0;

for i=1:numel(u_nom) % cost to go
    del_x=(x_des(i,:)-x_nom(i,:))';
    del_u=(u_des(i)-u_nom(i));
    cost=cost + sq(i,:) + del_x'*q(i,:)' + del_u'*r(i) + 0.5.*del_x'*Q*del_x + del_u'*R*del_u;
end
% cost;
J=h+cost; % total cost
% pause
end
