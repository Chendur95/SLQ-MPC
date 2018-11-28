function u_nom = SLQ_solve ( x_des,u_des,x_0,Q,R,Qf,t_vec,K, dt)
% Implementation of the SLQ alorithm from Fast nonlinear model predictive
% contro for unified trajectory optimization and tracking Neunert et al. ICRA 2016
%
%  Psuedocode --> page 1399 of the paper
% dimensionality check done


max_slq_iter=80;
max_lin_search_iter=100;
tol=0;iter=0;

slq_tol=0.05;
cost_tol=5;

x_des=x_des(t_vec,:);
u_des=u_des(t_vec,:);

n_steps=numel(t_vec);
% x_des=[x_des;[x_des(end,1)  x_des(end,2) ]];
% plot(x_des(:,1),x_des(:,2),'r.-');
% u_nom=0.01.*ones(n_steps,1);

for i=1:n_steps
    u_nom(i)=dot(K,(x_des(i,:)-x_0));
end

% x_0 supplied from MPC also equal to the first row of the x_nom
x_goal=[pi,0]; %% not needed

% P_tf, Q and R costs should come from MPC % we have 2 states so Q vectors are 2X2
    % only 1 control so R s are scalars
P_tf=Qf; p_tf=[10 10];
% these are matching dimansions from the paper equation
alpha=1;
alpha_d=2;             % linesearch update
u_max=2.0;

% system_physical_parameters= [1;1;9.8;0.3]; % [mass; length; accln. due to gravity; damping]

iter=0;
del_l=999;
del_u=999;
lt=zeros(n_steps,1); %  || del_l >=  slq_tol ||  iter< max_slq_iter  && 
while del_l >=  slq_tol  && iter< max_slq_iter
    l_old=lt;
    u_old=u_nom;
    x_nom=simulate_dynamics(u_nom, x_0,dt); %we do not need a x_nom tape as we generate it here by rollout
% to calculate updates for r, q  initialize
    sq=ones(n_steps,1); r=ones(n_steps,1);
        for i=1:n_steps
            sq(i,:)=diag(Q)'*(x_des(i,:)-x_nom(i,:))';
            r(i)=dot(diag(R),(u_des(i)-u_nom(i)));
        end
    % or change q and r to ones!
    q=ones(n_steps,2); %r=ones(n_steps,1);
    
    
    [Kt,lt]= riccati_like(P_tf, p_tf,Q,q,R,r,x_des,x_nom,t_vec);
    J_now=total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des,t_vec);
    %     pause
    J_new=inf;
    u_iter=0;  %  J_new -J_now >= cost_tol | && J_new  >= J_now ||
    while  u_iter<=max_lin_search_iter && J_new  > J_now 
        %     update control
        for i=numel(t_vec)
            u_nom(i)=u_old(i)+alpha*lt(i) + Kt(i,:)*(x_des(i,:)-x_nom(i,:))';
%             pause
            if u_nom(i)>=u_max
                u_nom(i)=u_max;
            elseif u_nom(i)<=-u_max
                u_nom(i)=-u_max;
            end
        end
        
        x_nom=simulate_dynamics(u_nom, x_0,dt);
        J_new=total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des,t_vec);
        alpha=alpha/alpha_d;
%         
        u_iter=u_iter+1;
    end
%     pause
    iter=iter+1;
% del_J=J_new-J_now    ;
    del_l=norm(l_old-lt,2);
% del_u=norm(u_old-u_nom,2);
    %     pause
end
%  plot to show how we did
% hold on
%  x_nom=simulate_dynamics(u_nom, x_0,dt);
%     plot(x_nom(:,1),x_nom(:,2),'b.-');
% plot(x_des(:,1),x_des(:,2),'g.-');
% figure(2)
% plot(u_nom)
% hold all
% pause
end


function [Kt_out,lt_out]= riccati_like(P_tf, p_tf,Q,q,R,r,x_des,x_nom,t_vec)
%  first populate the variable vectors
% Q  and q are the full state cost
% R and r are the control costs
% P is the quadratic difference cost
% p is the linear difference cost
n_steps=numel(t_vec); %2 x vector
P=zeros(n_steps,4); % we will reshape this to [2,2] before use
P(end,:)=[P_tf(1,1) P_tf(1,2) P_tf(2,1) P_tf(2,2)]; % from the input
p=zeros(n_steps,2); %will reshape this into [2,1]
p(end,:)=[p_tf(1) P_tf(2)];% from the input
m=1;l=1;g=1;b=1; % system properties

% initialize outputs
Kt_out=zeros(n_steps,2);
lt_out=zeros(n_steps,1);

for i = n_steps:-1:2
    Bt=[0;1]; %% 2X1 vector
    i;
    At= [0 1; -g/l*cos(x_nom(i,1)) -b]; %% 2X2 square matrix
    P_t_plus_1 = reshape(P(i,:),[2,2]); %% 2X2 square matrix
    p_t_plus_1 = p(i,:)';  %% 2X1 vector
    
    Ht=R+Bt'*P_t_plus_1 *Bt; %% scalar-- 1+ 1X2*2X2*2X1
    Gt=Bt'*P_t_plus_1*At ; %% 1X2 vector --1X2*2X2*2X2
    gt=r(i)+Bt'*p_t_plus_1;   %% scalar -- 1+1X2*2X1
    
    Kt=-inv(Ht)*Gt;  %% 1X2 vector--  1*1X2
    Kt_out(i,:)=Kt(:); %% for storage!-- we need to reshape them
    
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

Sim_x=zeros(n_steps,2); % as we have 2 states here
Sim_x(1,:)=x_0;
for i=1:n_steps-1
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
x_dot=[x(2);-g*sin(x(1))/l - b*x(2)/(m*l*l) + u/(m*l^2)];
end

function J= total_cost(x_nom,x_goal,u_nom,q,sq,Q,R,r,P_tf,p_tf,x_des,u_des,t_vec)
del_x_f=x_des(end,:)-x_nom(end,:);
% del_x_f=x_nom(end,:)-x_goal;
% x_nom(end,:)
h = del_x_f*p_tf' + 0.5.* del_x_f*P_tf*del_x_f'; %terminal stuff
cost=0;

for i=1:numel(t_vec) % cost to go
    del_x=(x_des(i,:)-x_nom(i,:))';
    del_u=(u_des(i)-u_nom(i));
    cost=cost + sq(i,:) + del_x'*q(i,:)' + del_u'*r(i) + 0.5.*del_x'*Q*del_x + del_u'*R*del_u;
end
J=h+cost; % total cost
end
