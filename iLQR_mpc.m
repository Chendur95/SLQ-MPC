clc
clear all
load('control_input.mat');
load('states_pendulum.mat'); 

x_des = store; 
u_des = u; 
 
R_lqr = 1;   
Q_lqr = 1500*eye(2); 
tot_steps = 200;
dt = 0.025;

plot(x_des(1,:), x_des(2,:), 'r.-'); 
title('Phase plot');
xlabel('Angle Theta');
ylabel('Angular velocity');
hold on; 
b=1; g=1; l=1; 
x_init = [0;0]; 
x_curr = x_init; 
B =[0;1]; 
n_steps = 5; 

% Look ahead 5 steps optimise and take first u 
control=[];  
flag = 1;
tic 
for t=1:1:(tot_steps-n_steps+2)
    xtracked(:,t) = x_curr; 
    if(t>(tot_steps - n_steps + 1))   % Correction for last n step steps
        n_steps = n_steps -1; 
        flag =2; 
    end
    t
    A = [0 1;
       -g*cos(x_curr(1))/l -b]; 
    K_lqr = lqr(A,B,Q_lqr,R_lqr);
    x_des_ilqr = x_des(: , (t+1) : (t + n_steps));
    u_des_ilqr = u_des(t:(t+n_steps-1));
    
    x_tf = x_des_ilqr(:,end); 
    A_tf = [0 1;
           -g*cos(x_tf(1))/l -b]; 
    [Kf,Qf] = lqr(A_tf,B, Q_lqr, R_lqr); 

    [u_ret] = ilqr_solve(x_des_ilqr, u_des_ilqr, x_curr, Qf, K_lqr, n_steps, dt, t);
    if(flag==1)
        u = u_ret(1);
        x = x_curr + dynamics(x_curr,u)'*dt ;
        x_curr = x;
    elseif(flag==2)
        kk = t; 
        steps = numel(u_ret);
        for k=1:1:steps
            x = x_curr + dynamics(x_curr,u_ret(k))'*dt ;
            x_curr = x;
            xtracked(:,kk) = x_curr; 
            kk = kk +1; 
        end
    end
    plot(xtracked(1,:), xtracked(2,:), 'b.-'); 
    legend('Desired trajectory', 'Tracked trajectory' ); 
   % pause;
end
toc
xtracked(:,tot_steps+1) = x_curr;
error = norm(xtracked - x_des)

function xdot = dynamics(x,u)
g = 1; l = 1; m = 1; b = 1;
xdot(:,1) = x(2); 
xdot(:,2) = -g*sin(x(1))/l - b*x(2)/(m*l*l) + u/(m*l*l); 
end 