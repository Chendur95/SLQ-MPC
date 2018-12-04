function u_ret = slq_solve_ImprovedSearch(x_des, u_des, x_st, Qf, K_lqr, Q, R, n_steps, dt)

P_tplus = Qf;                % P at time=10 
p_tplus = [10;10]; %P*x_nom(end,:)';
B = [0;1];
q = ones(n_steps,2)';    
r = ones(n_steps,1);      
u_max = 2;
u_min = -2;

del_l=999;
l=zeros(n_steps,1);

slq_tol = 0.05; 
iter =0;

u_nom = zeros(n_steps,1); 

for i=1:1:n_steps
   u_nom(i) = K_lqr*(x_des(:,i) - x_st);
end

while(iter<80 && del_l >= slq_tol)
    iter = iter +1;
    l_old = l; 
    x_curr = x_st; 
    x_nom(:,1) = x_st; 

    for i=2:(n_steps+1)       % If we have 100 u_nom values then till 101       
        xdot = dynamics(x_curr, u_nom(i-1))'; 
        x_nom(:,i) = x_curr + xdot*dt; 
        x_curr = x_nom(:,i);
    end

    for t = (n_steps):-1:1
    
        %u = u_nom(t); 
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
    x_new(:,1) = x_st; 
    x_curr = x_new(:,1); 
    J_new = Inf; 
    J_now = computeCost(pf,Pf,q,Q,R,r, x_nom, u_nom , x_des, u_des, n_steps); 
    uiter =0; 
    
    alpha = 1; 
    alpha_d = 1.4 ;
    flag = 1; 
    while(uiter<=100 && J_new>J_now && flag ~=0)
        uiter = uiter+1; 
        x_curr = x_st; 
        for rr =1:1:n_steps

            u_upd(rr) = u_nom(rr) + alpha*l(rr) + K(rr,:)*(x_des(:,rr)-x_nom(:,rr)); 
            
            if u_upd(rr)>=u_max
                u_upd(rr)=u_max;
            elseif u_upd(rr)<=u_min
                u_upd(rr)=u_min;
            end
            
            xdot = dynamics(x_curr, u_upd(rr))';
            x_new(:,rr+1) = x_curr + xdot*dt; 
            x_curr = x_new(:,rr+1); 
            a = 1; 
        end
        
        % J COMPUTE 
        Ptt = P(:,:,tt); %Send the last Pf value to find cost
        ptt = p(:,tt); 
        J_new = computeCost(ptt,Ptt,q,Q,R,r, x_new, u_upd, x_des, u_des, n_steps) ;
  %      alpha = alpha/alpha_d; 
  
          if(uiter==1)
            J_prev = J_new; 
          elseif (uiter==2)
            J_next = J_new; 
            if(J_next>J_prev)
               flag = 2;
            end
           J_prev = J_new; 
          else
            J_next = J_new; 
            if(J_next>=J_prev)
                flag = 0; 
            end
            J_prev = J_new; 
          end
       
          if(flag == 1)
            alpha = alpha/alpha_d;
          elseif(flag == 2)
            alpha = alpha*alpha_d; 
          end
  
    end
    del_J = J_new - J_now ;
    del_l=norm(l_old-l,2);
    u_nom = u_upd; 
end
u_ret = u_nom; 
end
%% Func to compute cost 

function Jcost = computeCost(p,P,q,Q,R,r,x_nom,u_nom, x_des, u_des, n_steps)

% Compute h
xf_del = -(x_nom(:,end) - x_des(:,end)); 
h_cost = xf_del'*p + 0.5*xf_del'*P*xf_del; 

% Compute L
L_cost = 0; 
for t=1:n_steps
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
