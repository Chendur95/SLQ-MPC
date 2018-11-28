clc
clear all
close all
load('control_input.mat');
load('states_pendulum.mat'); 
dt=0.025; 
n_steps=numel(u);

x_des= store;
u_des=u;
%hold on; 
x_init =[0;0]; 
u_nom = 2*u_des; 

Q = 50*eye(2); 
R = 1;  
Qf = 100*eye(2); 
P_tplus = Qf;                % P at time=10 
p_tplus = [10;10]; %P*x_nom(end,:)';
B = [0;1];
q = ones(200,2)';    
r = ones(200,1);      
u_max = 1.99; 

del_l=999;
l=zeros(n_steps,1);

slq_tol = 0.05; 
iter =0;
while(iter<80 && del_l >= slq_tol)
    iter = iter +1
    l_old = l; 
    x_curr =[0;0]; 
    x_nom(:,1) = [0;0]; 

    for i=2:(n_steps+1)       % If we have 100 u_nom values then till 101       
        xdot = dynamics(x_curr, u_nom(i-1))'; 
        x_nom(:,i) = x_curr + xdot*dt; 
        x_curr = x_nom(:,i);
    end

    plot(x_nom(1,:),x_nom(2,:),'b.-'); 
   % hold all; 
    for t = n_steps:-1:1
    
        u = u_nom(t); 
        A = findA(t,x_nom); 
    
        H = R + B'*P_tplus*B;   %scalar
        G = B'*P_tplus*A;       %1x2
    
        K(t,:) = -H\G;     %1x2
        g = r(t) + B'*p_tplus;   % Scalar The initialised p outisde loop is p(t+1)

        l(t) = -H\g - (H\G)*(x_des(:,t) - x_nom(:,t)); %scalar 
    
        p(:,t) = q(:,t) + A'*p_tplus + K(t,:)'*H*l(t) + K(t,:)'*g + G'*l(t); %2*1 vector - K*g - size ensure 
        P(:,:,t) = Q + A'*P_tplus*A + K(t,:)'*H*K(t,:) +K(t,:)'*G + G'*K(t,:); 
        
        p_tplus = p(:,t);   %2x1
        P_tplus = P(:,:,t); %2x2 
    end
    
    tt = n_steps +1; 
    P(:,:,tt) = Qf;         
    p(:,tt) = [1;1];       
    pf = p(:,tt);
    Pf = P(:,:,tt); 
    
    % LINEAR SEARCH - U_NOM - UPDATED 
    x_new(:,1) = x_init; 
    x_curr = x_new(:,1); 
    J_new = Inf; 
    J_now = computeCost(pf,Pf,q,Q,R,r, x_nom, u_nom , x_des, u_des); 
    uiter =0; 
    
    alpha = 2; 
    alpha_d = 10 ;
    %del_J = -10;  %Just to enter the loop 
    
    while(uiter<100 && J_new>J_now)
        uiter = uiter+1; 
        x_curr = x_init; 
        for rr =1:1:n_steps

            u_upd(rr) = u_nom(rr) + alpha*l(rr) + K(rr,:)*(x_des(:,rr)-x_nom(:,rr)); 
            
            if u_upd(rr)>=u_max
                u_upd(rr)=u_max;
            elseif u_upd(rr)<=-u_max
                u_upd(rr)=-u_max;
            end
            
            xdot = dynamics(x_curr, u_upd(rr))';

            x_new(:,rr+1) = x_curr + xdot*dt; 
            x_curr = x_new(:,rr+1); 
            a = 1; 
        end
        
        % J COMPUTE 
        Ptt = P(:,:,tt); %Send the last Pf value to find cost
        ptt = p(:,tt); 
        J_new = computeCost(ptt,Ptt,q,Q,R,r, x_new, u_upd, x_des, u_des) ;
        alpha = alpha/alpha_d; 
      %  u_nom = u_upd; 
    end
    J_new
    del_J = J_new - J_now 
    del_l=norm(l_old-l,2)
    u_nom = u_upd; 
    pause
end
hold on;
    plot(x_des(1,:),x_des(2,:),'r.-');

%% Func to compute cost 

function Jcost = computeCost(p,P,q,Q,R,r,x_nom,u_nom, x_des, u_des)

% Compute h
xf_del = -(x_nom(:,end) - x_des(:,end)); 
h_cost = xf_del'*p + 0.5*xf_del'*P*xf_del; 

% Compute L
L_cost = 0; 
for t=1:200
    u_del = (u_des(t) - u_nom(t)); 
    x_del = (x_des(:,t) - x_nom(:,t));  
    L_cost = L_cost + x_del'*q(:,t) +u_del'*r(t) + 0.5*x_del'*Q*x_del +u_del'*R*u_del; 
end

Jcost = h_cost + L_cost; 
end

%% Func to Linearize dynamics and find A matrix  

function A = findA(t,x_nom)    % t is current Time step 
theta = x_nom(1,t); 
b = 1; g = 1; l=1; 
A = [0 1;
    -g*cos(theta)/l -b]; 
end

function xdot = dynamics(x,u)
g = 1; l = 1; m = 1; b = 1;
xdot(:,1) = x(2); 
xdot(:,2) = -g*sin(x(1))/l - b*x(2)/(m*l*l) + u/(m*l*l); 

end
