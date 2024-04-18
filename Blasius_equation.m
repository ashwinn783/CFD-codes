clc, clear, close
tic
%solving the blasius equation
%parameters
unk=0.0005;%delta_eta
length_eta=8;
n_eta=length_eta/unk;
delta_x=0.0005;
length_x=1;
nx=length_x/delta_x;
nu=1e-3;
u_infinity=10;

%numerical solution of blasius equation
[guess]=Iteration_initial_condition_Blasius_equation;
[y_1,y_2,y_3]=Runge_kutta_method(guess,unk);

%construct matrix
eta=0:0.0005:8;
y_00005=zeros(length(eta),1);
u_00005=zeros(length(eta),1);
v_00005=zeros(length(eta),1);
y_05=zeros(length(eta),1);
u_05=zeros(length(eta),1);
v_05=zeros(length(eta),1);

%horizontal and vertical velocity
x=0.0005;
y_00005(:,1)=eta./(length_x*(u_infinity/(nu*x))^0.5);
u_00005(:,1)=y_2;
v_00005(:,1)=0.5*(nu/(u_infinity*x))^0.5.*(eta.*y_2-y_1);
x=0.5;
y_05(:,1)=eta./(length_x*(u_infinity/(nu*x))^0.5);
u_05(:,1)=y_2;
v_05(:,1)=0.5*(nu/(u_infinity*x))^0.5.*(eta.*y_2-y_1);

figure('Name','velocity profile at x=0.0005 and 0.5','NumberTitle','off')

hold on

title('Blasius equation')
plot(y_00005,u_00005,'b--','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k')
xlim([0 0.1])
ylim([0 1]);

plot(y_05,u_05,'b-','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k')
xlim([0 0.1])
ylim([0 1]);

plot(y_00005,v_00005,'r--','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k')
xlim([0 0.1])
ylim([0 1]);

plot(y_05,v_05,'r-','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k')
xlim([0 0.1])
ylim([0 1]);

Legend1={'u profile at x=0.0005','u profile at x=0.5','v profile at x=0.0005','v profile at x=0.5'};
legend(Legend1,'FontSize',15);
grid on;

hold off

toc


