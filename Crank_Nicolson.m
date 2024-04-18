clc, clear, close
tic
%Solving the 2D heat equation using the Crank-Nicolson scheme
%% Calculation
%Parameters
delta_x=0.025;
length_x=1;
nx=length_x/delta_x; %Number of grid intervals along the x-direction
x=0:delta_x:1;
delta_y=0.025;
length_y=1;
ny=length_y/delta_y; %Number of grid intervals along the y-direction
y=0:delta_y:1;
delta_t=0.0001;
range_t=0.16;
nt=range_t/delta_t; %Number of time intervals
t=0:delta_t:0.16;

%Coefficients in the equation
alpha=delta_t/(2*delta_x^2);
fac_front=1+2*alpha;
fac_back=1-2*alpha;
fac_1bc=0.5*(1-2*alpha);
fac_2bc=0.5*(1+2*alpha);

%Assigning uniform grid intervals and constructing coefficient matrix
first_row=eye(1,nx+1);
last_row=[zeros(1,nx) 1];
tridiag=(fac_front*eye(length(x))-alpha*diag(ones(1,length(x)-1),-1)-alpha*diag(ones(1,length(x)-1),1));
tridiag(1,:)=first_row;
tridiag(end,:)=last_row;
coe_matrix=tridiag;

%Constructing the initial temperature matrix
for i=1:nx+1
    for j=1:ny+1
        w(i,j)=0;
    end
end

%Construct the asterisk matrix
for i=1:nx+1
    for j=1:ny+1
        w_asterisk(i,j)=0;
    end
end

%Initial and boundary conditions
top_bc=1-y'.^3;
bottom_bc=1-sin(pi/2*y');
right_bc=0;
left_bc=1;
w_0=zeros(nx+1,1);
time_t0=0;
c=nx+1;
B=zeros(nx+1,1);

%Inserting boundary conditions into the temperature matrix
w(:,end)=right_bc;
w(:,1)=left_bc;
w(end,:)=bottom_bc;
w(1,:)=top_bc;

%Constructing the stored matrix for time evolution
time_ev=zeros(1,nt+1);
temp_ev=zeros(1,nt+1);
time_ev(1,1)=0;
temp_ev(1,1)=0;
time=time_t0;

%Crank-Nicolson scheme
[X,Y]=meshgrid(x,y);
%Transfer matrix notation 
for time_step=1:nt
    for j=1:ny+1
        i=1:nx+1;
        if time_step==1
            B=w_0;
            R=coe_matrix\B;
            w_asterisk(:,j)=R;
        elseif j==1
            B(1,1)=fac_1bc*w(1,j)+fac_2bc*w(1,j)+0.5*alpha*(w(1,j+1)+-w(1,j+1));
            B(nx+1,1)=fac_1bc*w(nx+1,j)+fac_2bc*w(nx+1,j)+0.5*alpha*(w(nx+1,j+1)-w(nx+1,j+1));
            B(i,1)=fac_back*w(i,j)+alpha*w(i,j+1);
            R=coe_matrix\B;
            w_asterisk(:,j)=R;
        elseif j==ny+1
            B(1,1)=fac_1bc*w(1,j)+fac_2bc*w(1,j)+0.5*alpha*(w(1,j-1)-w(1,j-1));
            B(nx+1,1)=fac_1bc*w(nx+1,j)+fac_2bc*w(nx+1,j)+0.5*alpha*(w(nx+1,j-1)-w(nx+1,j-1));
            B(i,1)=fac_back*w(i,j)+alpha*w(i,j-1);
            R=coe_matrix\B;
            w_asterisk(:,j)=R;
        else 
            B(1,1)=fac_1bc*w(1,j)+fac_2bc*w(1,j)+0.5*alpha*(w(1,j+1)+w(1,j-1)-w(1,j+1)-w(1,j-1));
            B(nx+1,1)=fac_1bc*w(nx+1,j)+fac_2bc*w(nx+1,j)+0.5*alpha*(w(nx+1,j+1)+w(nx+1,j-1)-w(nx+1,j+1)-w(nx+1,j-1));
            B(i,1)=fac_back*w(i,j)+alpha*w(i,j+1)+alpha*w(i,j-1);
            R=coe_matrix\B;
            w_asterisk(:,j)=R;
        end
    end
    for i=2:nx
        j=2:ny;
        B(1,1)=1;
        B(nx+1,1)=0;
        B(j,1)=fac_back*w_asterisk(i,j)+alpha*w_asterisk(i+1,j)+alpha*w(i-1,j);
        R=coe_matrix\B;
        w(i,:)=R.';
    end
    w(:,end)=right_bc;
    w(:,1)=left_bc;
    w(end,:)=bottom_bc;
    w(1,:)=top_bc;
    temp_ev(1,time_step+1)=w(17,17);
    time=time+delta_t;
    time_ev(1,time_step+1)=time;
  %{
  if time_step==100
        figure('Name','Entire domain t=0.01','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    elseif time_step==200
        figure('Name','Entire domain t=0.02','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    elseif time_step==400
        figure('Name','Entire domain t=0.04','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    elseif time_step==800
        figure('Name','Entire domain t=0.08','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    elseif time_step==1600
        figure('Name','Entire domain t=0.16','NumberTitle','off')
        surf(x,y,w)
        axis([0 1 0 1 0 1])
        xlabel('x')
        ylabel('y')
        zlabel('T');
    end
  %}
end

%% Plotting time evolution of temperature at x = y = 0.4
figure('Name','Time evolution','NumberTitle','off')
plot(time_ev,temp_ev,'k-','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('Time','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('Temperature','Fontsize',15,'FontWeight','bold','Color','k');
grid on;

%% Plotting vertical temperature profile at t = 0.16 and x = 0.4
figure('Name','Vertical temperature','NumberTitle','off')
plot(y,w(17,:),'k-','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('Temperature','Fontsize',15,'FontWeight','bold','Color','k');
grid on;
toc
