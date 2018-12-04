function u_ret = ilqr_solve(x_des, u_des, x_st, Qf, K_lqr, num_steps, dt, t)

Q = 100*eye(2); 
R=1;  
x_init = x_st;  

B = [0;1];  % B matrix of linearised dynamics 
tspan =linspace(dt*(num_steps+t), dt*t, num_steps);     
S2f = Qf;
times = linspace(dt*t, dt*(num_steps+ t), num_steps);       

for i=1:1:num_steps
   u_nom(i) = K_lqr*(x_des(:,i) - x_st);
end

x_curr = x_st; 

for i=2:(num_steps+1)       % If we have 100 u_nom values then till 101       
    xdot = dynamics_rollout(x_curr, u_nom(i-1))'; 
    x_nom(:,i-1) = x_curr + xdot*dt; 
    x_curr = x_nom(:,i-1);
end

x_des = x_des';
x_nom = x_nom' ;

ilqr_iter = 3; 
% iLQR 
u=[];
for i=1:1:ilqr_iter
    [t, S2] = ode45(@(t,S2) riccatiS2(t,S2,x_nom, times,B,Q),tspan,S2f); 
    tS2 = t; 
    S1f = -2*Qf*(x_des(end,:)'-x_nom(end,:)'); 

    [t2, S1] = ode45(@(t2,S1) riccatiS1(t2,S1,x_nom,u_nom,x_des,u_des, times,B,Q,R,S2,tS2),tspan,S1f); 
    
    tS2 = flip(t); 
    S2 = flipud(S2); 
    tS1 = flip(t2);
    S1 = flipud(S1); 

    for j = 1:num_steps
        s2=reshape(S2(j,:),2,2);
        s1=reshape(S1(j,:),2,1);
        M2(j,:) = B'*s2;
        M1(j,:) = B'*0.5*s1;
    end
    
    [t3, x] = ode45(@(t3,x) dynamics(t3, x, x_nom, x_des, u_nom, u_des, times, B,Q, M1, M2, tS1, tS2), times, x_init);  
    
    for k=1:num_steps
        u(k) = u_des(k) - M2(k,:)*(x(k,:)-x_nom(k,:))' - M1(k);
    end
     x_nom = x; 
     u_nom = u; 
    Error = norm(x_des-x_nom);
end
u_ret = u_nom; 
end
%% Riccati functions 

function dS_two_dt = riccatiS2(t,S2,x_nom, T,B,Q) 
A = findA(t,T,x_nom); 
S2 = reshape(S2,size(Q));
dS_two_dt = -(Q-S2*B*transpose(B)*S2 +S2*A+transpose(A)*S2); 
dS_two_dt = dS_two_dt(:); 
end

function dS_one_dt = riccatiS1(t,S1,x_nom, u_nom, x_des, u_des, T,B,Q,R,S2, tS2) 
A = findA(t,T,x_nom);  
S1 = reshape(S1,size(B)); 
ind1 = interp1(tS2, S2, t);  
S2 = reshape(ind1, size(Q)); 
x_d = interp1(T, x_des, t); 
x_n = interp1(T, x_nom, t);
u_d = interp1(T, u_des, t); 
u_n = interp1(T, u_nom, t);
udash = u_d - u_n; 
xdash = (x_d - x_n)'; 
dS_one_dt = -(-2*Q*xdash +(transpose(A)-S2*B*transpose(B))*S1+2*S2*B*udash);  
dS_one_dt = dS_one_dt(:); 
end

%% Dynamics function 
function xdot = dynamics(t, x, x_nom, x_des,u_nom, u_des, T, B,Q, M1, M2, tS1, tS2) 
b = 1; 
M2_curr = interp1(tS2, M2, t);  
M1_curr = interp1(tS1, M1, t); 

u_d = interp1(T, u_des, t);
u_n = interp1(T, u_nom, t); 
x_d = interp1(T, x_des, t); 
x_n = interp1(T, x_nom, t); 
xdash = x - x_n'; 
ustar = u_d - M2_curr*xdash - M1_curr; 
xdot = [x(2); ustar-b*x(2)-sin(x(1))]; 
end

%% Func to Linearize dynamics and find A matrix  
function A = findA(t,T,x_nom)
theta =interp1(T, x_nom(:,1), t); 
b = 1; 
A = [0 1;
    -cos(theta) -b]; 
end

function xdot = dynamics_rollout(x,u)
g = 1; l = 1; m = 1; b = 1;
xdot(:,1) = x(2); 
xdot(:,2) = -g*sin(x(1))/l - b*x(2)/(m*l*l) + u/(m*l*l); 
end

