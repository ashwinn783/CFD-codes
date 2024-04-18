clc, clear, close
tic
%solving the boundary layer equation
%% calculation
%parameters
delta_x=0.0005;%grid distance along x axis
length_x=1;
nx=length_x/delta_x;
x=0:delta_x:1;
delta_y=0.0005;%grid diatance along y axis
Re=1e4;
delta=5/(Re)^0.5;
length_y=2*delta;
ny=length_y/delta_y;
y=0:delta_y:length_y;
u_infinity=1;

%construct velocity matrix
u=zeros(ny+1,nx+1);
v=zeros(ny+1,nx+1);

%boundary conditions
u(:,1)=1;%left_bc
v(:,1)=0;
u(ny+1,:)=0;%bottom_bc
v(ny+1,:)=0;
u(1,:)=u_infinity;%top_bc
v(1,:)=0;

%coefficient of discretization equation
co_1=delta_x/delta_y^2;
co_2=delta_x/(2*delta_y);
co_3=delta_y/(2*delta_x);
co_4=1/Re;


%numerical method to solve boundary layer equation
[X,Y]=meshgrid(x,y);
for i=1:nx
     for j=2:ny
         u(j,i+1)=u(j,i)+co_4*co_1*(u(j+1,i)-2*u(j,i)+u(j-1,i))*(1/u(j,i))-co_2*(v(j,i)/u(j,i))*(u(j+1,i)-u(j-1,i));
         v(j,i+1)=v(j-1,i+1)-co_3*(u(j,i+1)-u(j,i)+u(j,i+1)-u(j-1,i));
     end
end
toc

%% figure
% the whole domain
figure('Name','u of the whole domain','NumberTitle','off')
y_reality=fliplr(y);
contourf(x,y_reality,u)
set(gca,'YScale','log')
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('x','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('y','Fontsize',15,'FontWeight','bold','Color','k')
a=colorbar;
ylabel(a,'u','FontSize',15)
caxis([0 1]);


figure('Name','v of the whole domain','NumberTitle','off')
y_reality=fliplr(y);
contourf(x,y_reality,v)
set(gca,'YScale','log')
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('x','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('y','Fontsize',15,'FontWeight','bold','Color','k')
a=colorbar;
ylabel(a,'v','FontSize',15)
caxis([0 0.2]);

% velocity profile
figure('Name','velocity profile at x=0.0005 and 0.5','NumberTitle','off')
hold on

title('Boundary layer equation')
plot(y_reality,u(:,2),'b--','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k');

plot(y_reality,u(:,0.5/delta_x+1),'b-','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k')

plot(y_reality,v(:,2),'r--','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k');

plot(y_reality,v(:,0.5/delta_x+1),'r-','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k');

Legend={'u profile at x=0.0005','u profile at x=0.5','v profile at x=0.0005','v profile at x=0.5'};
legend(Legend,'FontSize',15);
grid on;

hold off



