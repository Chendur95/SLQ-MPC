function plot_cart_pole (State,Time)
pos = State(1,:);
angle = pi-State(2,:);

%Position of the cart over time:
Cart = [pos;zeros(size(pos))];

%Position of the pendulum bob over time:  (Assume that length==1)
Bob = Cart + [-sin(angle);cos(angle)];

%Pack up for plotting:
x = [Cart;Bob];
t = Time;

%Figure out the extents of the axis
horizontalAll = [Cart(1,:), Bob(1,:)];
verticalAll = [Cart(2,:), Bob(2,:)];
param.axis = [min(horizontalAll),max(horizontalAll),...
  min(verticalAll),max(verticalAll)];
param.clearFigure = true;

%Pass the trace of the Bob's path to the plotting function.
param.Bob = Bob;

%Set up parameters for animation:
P.plotFunc = @(t,x)plotPendulumCart(t,x,param);
P.speed = 1.5;
P.figNum = gcf;

%Call the animation function
Animate(t,x,P)
end

function plotPendulumCart(t,x,param)

Cart_Width = 0.2;
Cart_Height = 0.1;
Limits = param.axis;

Cart = x(1:2);
Bob = x(3:4);

if param.clearFigure    %THEN RUNNING ANIMATION
    clf;
    hold on;
    
    %Plot Cart
    x = Cart(1) - 0.5*Cart_Width;
    y = -0.5*Cart_Height;
    w = Cart_Width;
    h = Cart_Height;
    h = rectangle('Position',[x,y,w,h],'LineWidth',4,'Curvature',[0.3,0.3]);
    set(h,'EdgeColor',[0.1,0.8,0.1])
    
    %Plot Pendulum
    Rod_X = [Cart(1), Bob(1)];
    Rod_Y = [Cart(2), Bob(2)];
    plot(Rod_X,Rod_Y,'k-','LineWidth',4)
    
    %Plot Bob
    plot(Bob(1),Bob(2),'k.','MarkerSize',50)
    
    %Plot Rails
    plot([Limits(1) Limits(2)],-0.5*Cart_Height*[1,1],'k-','LineWidth',2)
    
    %Title
    title(['Simulation Time: ' num2str(t) ' s'])
    
    
    %These commands keep the window from automatically rescaling in funny ways.
    axis(Limits);
    axis('equal');
    axis manual;
    axis off;
    
    
else    %THEN RUNNING STOP ACTION
    
    hold on;
    
    %Plot Cart
    x = Cart(1) - 0.5*Cart_Width;
    y = -0.5*Cart_Height;
    w = Cart_Width;
    h = Cart_Height;
    rectangle('Position',[x,y,w,h],'LineWidth',2);
    
    %Plot Pendulum
    Rod_X = [Cart(1), Bob(1)];
    Rod_Y = [Cart(2), Bob(2)];
    plot(Rod_X,Rod_Y,'k-','LineWidth',2)
    
    %Plot Bob and hinge
    plot(Bob(1),Bob(2),'k.','MarkerSize',30)
    plot(Cart(1),Cart(2),'k.','MarkerSize',30)
    
    %Plot Rails
    plot([Limits(1) Limits(2)],-0.5*Cart_Height*[1,1],'k-','LineWidth',2)
    
    %Plot trace of the Bob's path:
    plot(param.Bob(1,:),param.Bob(2,:),'k:','LineWidth',1)
    
    
    %These commands keep the window from automatically rescaling in funny ways.
    axis(Limits);
    axis('equal');
    axis manual;
    axis off;
    
end
end
 
function Animate(t,x,P)
%Animate(t,x,P)
%
%FUNCTION:
%   Animate is used to animate a system with state x at the times in t.
%
%INPUTS:
%   t = [1xM] vector of times, Must be monotonic: t(k) < t(k+1)
%   x = [NxM] matrix of states, corresponding to times in t
%   P = animation parameter struct, with field:
%     plotFunc = @(t,x) = function handle to create a plot
%       t = a scalar time
%       x = [Nx1] state vector
%     speed = scalar multiple of time, for playback speed 
%
%OUTPUTS: 
%   Animation based on data in t and x.
%

if ~isfield(P,'figNum')
    P.figNum=1000;  %Default to figure 1000
end
figure(P.figNum)


Loop_Time = 0;    %store how long has the simulation been running
T_end = t(end);   %Ending time of one step

t = t-min(t);  %Force things to start at t=0;

tic;    %Start a timer
while Loop_Time < T_end;  %Loop while the CPU time is less than the end of the simulation's time
    
    %Interpolate to get the new point:
    size(t)
    size(x)
    xNow = interp1(t',x',Loop_Time,'pchip','extrap')';
    
    %Call the plot command
   feval(P.plotFunc,Loop_Time,xNow);
   drawnow;
   
       %The next few lines pull the time step that is closest to the real time
    RealTime = toc;
    Loop_Time = RealTime*P.speed;   %Get the current time  (Taking slow motion into accunt if desired)

end

end