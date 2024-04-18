clc, clear, close
tic
%Solving the 2D heat equation with explicit Euler scheme
%% Calculation
%Parameters
delta_x=0.025;
length_x=1;
nx=length_x/delta_x; %Number of grid intervals along x-direction
x=0:delta_x:1;
delta_y=0.025;
length_y=1;
ny=length_y/delta_y; %Number of grid intervals along y-direction
y=0:delta_y:1;
delta_t=0.0001;
range_t=0.16;
nt=range_t/delta_t; %Number of time intervals
t=0:delta_t:0.16;

%Constructing coefficient matrces Ax and Ay
Ax=(-2*eye(length(x))+diag(ones(1,length(x)-1),-1)+diag(ones(1,length(x)-1),1));
Ay=(-2*eye(length(y))+diag(ones(1,length(y)-1),-1)+diag(ones(1,length(y)-1),1));

%Initial and boundary conditions
top_bc=1-y'.^3;
bottom_bc=1-sin(pi/2*y');
right_bc=0;
left_bc=1;
w_t0=0;
time_t0=0;

%Constructing the stored matrix for time evolution
time_ev=zeros(1,nt+1);
temp_ev=zeros(1,nt+1);
time_ev(1,1)=0;
temp_ev(1,1)=0;

%Expilcit Euler scheme
[X,Y]=meshgrid(x,y);
w=w_t0;
time=time_t0;
alpha=1/(delta_x)^2; %Uniform delta_x and delta_y
for n=1:nt
    w=w+(alpha*w*Ax+alpha*Ay*w)*delta_t; %Explicit Euler equation
    w(:,ny+1)=right_bc;
    w(:,1)=left_bc;
    w(end,:)=bottom_bc;
    w(1,:)=top_bc;
    temp_ev(1,n+1)=w(17,17);
    time=time+delta_t;
    time_ev(1,n+1)=time;
    if n==100
        figure('Name','Entire domain t=0.01','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    elseif n==200
        figure('Name','Entire domain t=0.02','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    elseif n==400
        figure('Name','Entire domain t=0.04','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    elseif n==800
        figure('Name','Entire domain t=0.08','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    elseif n==1600
        figure('Name','Entire domain t=0.16','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    end
end

%% Plotting time evolution of temperature at x=y=0.4
figure('Name','Time evolution','NumberTitle','off')
plot(time_ev,temp_ev,'k-','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('Time','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('Temperature','Fontsize',15,'FontWeight','bold','Color','k');
grid on;

%% Plotting vertical temperature profile at t=0.16 and x=0.4
figure('Name','Vertical temperature','NumberTitle','off')
plot(y,w(17,:),'k-','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('Temperature','Fontsize',15,'FontWeight','bold','Color','k');
grid on;
toc
