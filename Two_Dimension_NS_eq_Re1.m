%% Reference: https://www.youtube.com/watch?v=iyqAib3QXY8&t=869s
clc, clear, close
tic
%solving the 2D Navier stokes equation for lid driven flow
%% Calculation
%parameters
N = 51; %the number of nodes 
length_x = 1;
h = length_x/(N-1); %equidistant:h=dx=dy
delta_t = 1e-3; 
t=0; % time
x = 0:h:length_x; %X domain span
y = 0:h:length_x; %Y domain span
Re = 1; %Reynolds number
nu = 1/Re;
%correction factor
alpha = 0.6; % for velocity
alpha_p = 0.3; % for pressure
%% Initializing the variables

%Final collocated variables
u_final(N,N)=0;
v_final(N,N)=0;
p_final(N,N)=1;
u_final(1,:) = 1;

%Staggered variables
u(N+1,N)=0; % u velocity grid has an additional row, add a row to u
u_asterisk(N+1,N)=0;
d_e(N+1,N)=0;
a_e_u(N,N)=-nu;
a_w_u(N,N)=-nu;
a_n_u(N,N)=-nu;
a_s_u(N,N)=-nu;
a_ij_u(1:N,1:N)=h^2/delta_t+4*nu;
v(N,N+1)=0; % v velocity has an additional column, add a row to v
v_asterisk(N,N+1)=0;
d_n(N,N+1)=0;
a_e_v(N,N)=-nu;
a_w_v(N,N)=-nu;
a_n_v(N,N)=-nu;
a_s_v(N,N)=-nu;
a_ij_v(1:N,1:N)=h^2/delta_t+4*nu;
p(N+1,N+1)=0; % pressure grid is one point larger than the physical grid on each side
p_star(N+1,N+1)=0;
p_prime(N+1,N+1)=0; % Pressure correction p'
b(N+1,N+1)=0;
a_e_p(N,N)=0;
a_w_p(N,N)=0;
a_n_p(N,N)=0;
a_s_p(N,N)=0;
a_ij_p(1:N,1:N)=0;

u(1,:)=2; % defining u(1,:)=1 and u(2,:)=1 has the same effect

u_corrected(N+1,N)=0;
v_corrected(N,N+1)=0;
p_corrected(N+1,N+1)=1;
u_corrected(1,:)=2;

% Cells to store variable values at given time point
cell_u=cell(1,6);
cell_v=cell(1,6);
cell_p=cell(1,6);
sum=0;

%% Solving the governing equations
error = 1;
iterations = 0;
error_req = 1e-7; %final required error residual


while error > error_req
    % x-momentum eq. - Interior
    for j = 2:N
        for i = 2:N - 1
            u_e = 0.5*(u(j,i) + u(j,i+1));  % known velocities
            u_w = 0.5*(u(j,i) + u(j,i-1));
            v_n = 0.5*(v(j-1,i) + v(j-1,i+1));
            v_s = 0.5*(v(j,i) + v(j,i+1));
            
            a_e_u(j,i) = 0.5*u_e*h - nu;
            a_w_u(j,i) = -0.5*u_w*h - nu;
            a_n_u(j,i) = 0.5*v_n*h - nu;
            a_s_u(j,i) = -0.5*v_s*h - nu;
            
            a_ij_u(j,i) = h^2/delta_t + 0.5*u_e*h - 0.5*u_w*h + 0.5*v_n*h - 0.5*v_s*h + 4*nu;
            
            A_e = -h;
            d_e(j,i) = A_e/a_ij_u(j,i);
            
            u_asterisk(j,i) = -(a_e_u(j,i)*u(j,i+1) + a_w_u(j,i)*u(j,i-1) + a_n_u(j,i)*u(j-1,i) + a_s_u(j,i)*u(j+1,i))/a_ij_u(j,i) + d_e(j,i)*(p(j,i+1) - p(j,i)) + u(j,i)*h^2/(delta_t*a_ij_u(j,i));
        end
    end
    
    % x-momentum eq. - Boundary
    u_asterisk(1,:) = 2 - u_asterisk(2,:); % top veclocity of 1
    u_asterisk(N + 1,:) = -u_asterisk(N,:); %
    u_asterisk(2:N,1) = 0; % left wall has u velocity 0
    u_asterisk(2:N,N) = 0; % right wall has u velocity 0
    
    % y-momentum eq. - Interior
    for j = 2:N - 1
        for i = 2:N
            u_e = 0.5*(u(j,i) + u(j+1,i));
            u_w = 0.5*(u(j,i-1) + u(j+1,i-1));
            v_n = 0.5*(v(j-1,i) + v(j,i));
            v_s = 0.5*(v(j,i) + v(j+1,i));
            
            a_e_v(j,i) = 0.5*u_e*h - nu;
            a_w_v(j,i) = -0.5*u_w*h - nu;
            a_n_v(j,i) = 0.5*v_n*h - nu;
            a_s_v(j,i) = -0.5*v_s*h - nu;
            
            a_ij_v(j,i) = h^2/delta_t + 0.5*u_e*h - 0.5*u_w*h + 0.5*v_n*h - 0.5*v_s*h + 4*nu;
            
            A_n = -h;
            d_n(j,i) = A_n/a_ij_v(j,i);
            
            v_asterisk(j,i) = -(a_e_v(j,i)*v(j,i+1) + a_w_v(j,i)*v(j,i-1) + a_n_v(j,i)*v(j-1,i) + a_s_v(j,i)*v(j+1,i))/a_ij_v(j,i) + d_n(j,i)*(p(j,i) - p(j+1,i)) + v(j,i)*h^2/(delta_t*a_ij_v(j,i));
        end
    end
    
    % y-momentum eq. - Boundary
    v_asterisk(:,1) = -v_asterisk(:,2);
    v_asterisk(:,N + 1) = -v_asterisk(:,N);
    v_asterisk(1,2:N) = 0;
    v_asterisk(N,2:N) = 0;
    
    % Zeroing the corrections to begin with (for each loop; if values converges, pc should be zero)
    p_prime(1:N+1,1:N+1)=0;
    
    % Continuity equation a.k.a. pressure correction - Interior
    for j = 2:N
        for i = 2:N
            a_e_p(j,i) = -h^2/(h^2/delta_t+(a_ij_u(j,i)-h^2/delta_t));
            a_w_p(j,i) = -h^2/(h^2/delta_t+(a_ij_u(j,i-1)-h^2/delta_t));
            a_n_p(j,i) = -h^2/(h^2/delta_t+(a_ij_v(j-1,i)-h^2/delta_t));
            a_s_p(j,i) = -h^2/(h^2/delta_t+(a_ij_v(j,i)-h^2/delta_t));
            a_ij_p(j,i) = -(a_e_p(j,i) + a_w_p(j,i) + a_n_p(j,i) + a_s_p(j,i));
            b(j,i) = (u_asterisk(j,i) - u_asterisk(j,i-1))*h - (v_asterisk(j,i) - v_asterisk(j-1,i))*h;
            
            p_prime(j,i) = -(a_e_p(j,i)*p_prime(j,i+1) + a_e_p(j,i)*p_prime(j,i-1) + a_e_p(j,i)*p_prime(j-1,i) + a_e_p(j,i)*p_prime(j+1,i) + b(j,i))/a_ij_p(j,i);
        end
    end
    
    % Correcting the pressure field
    for j = 2:N
        for i = 2:N
            p_corrected(j,i) = p(j,i) + alpha_p*p_prime(j,i);
        end
    end
    
    % Continuity eq. - Boundary (no permeation, pressure across the wall should be equal to each other)
    p_corrected(1,:) = p_corrected(2,:);
    p_corrected(N + 1,:) = p_corrected(N,:);
    p_corrected(:,1) = p_corrected(:,2);
    p_corrected(:,N + 1) = p_corrected(:,N);
    
    % Correcting the velocities
    for j = 2:N
        for i = 2:N - 1
            u_corrected(j,i) = u_asterisk(j,i) + alpha*((-h/a_ij_u(j,i))*(p_prime(j,i+1) - p_prime(j,i)));
        end
    end
    
    % x-momentum eq. - Boundary
    u_corrected(1,:) = 2 - u_corrected(2,:);
    u_corrected(N + 1,:) = -u_corrected(N,:);
    u_corrected(2:N,1) = 0;
    u_corrected(2:N,N) = 0;
    
    for j = 2:N - 1
        for i = 2:N
            v_corrected(j,i) = v_asterisk(j,i) + alpha*(-h/a_ij_v(j,i))*(p_prime(j,i) - p_prime(j+1,i));
        end
    end
    
    % y-momentum eq. - Boundary
    v_corrected(:,1) = -v_corrected(:,2);
    v_corrected(:,N + 1) = -v_corrected(:,N);
    v_corrected(1,2:N) = 0;
    v_corrected(N,2:N) = 0;
    
    
    % Continuity residual as error measure
    error = 0;
    for j = 2:N
        for i = 2:N
            error = error + abs(b(j,i));
        end
    end
    
    % After the converged solution, we map the staggered variables to
    % collocated variables
    for j = 1:N
        for i = 1:N
            u_final(j,i) = 0.5*(u(j,i) + u(j+1,i));
            v_final(j,i) = 0.5*(v(j,i) + v(j,i+1));
            p_final(j,i) = 0.25*(p(j,i) + p(j,i+1) + p(j+1,i) + p(j+1,i+1));
        end
    end
    
    % Generate figures for t=0:0.02:0.1
    if t<=0.1+1e-5 && t>=2e-2
        if (rem(t,0.02)) <=1e-6
            sum=sum+1;
            cell_u(1,sum)={u_final};
            cell_v(1,sum)={v_final};
            cell_p(1,sum)={p_final};
        end
    end
    
    % Finishing the iteration
    u = u_corrected;
    v = v_corrected;
    p = p_corrected;
    iterations = iterations + 1;
    t=t+delta_t;
end

%% figure
x_dom = ((1:N)-1).*h;
y_dom = 1-((1:N)-1).*h;
[X,Y] = meshgrid(x_dom,y_dom);
u_05=zeros(51,1);
u_05=u_final(:,26);
v_05=zeros(51,1);
v_05=v_final(:,26);
y2=flip(y);



figure(2);
contourf(X,Y,v_final, 21, 'LineStyle', 'none')
colorbar
a=colorbar;
ylabel(a,'v')
colormap('jet')
title('v velocity of whole domain at Re=0.1, dt=0.001')
xlabel('x')
ylabel('y')

figure(3);
contourf(X,Y,u_final, 21, 'LineStyle', 'none')
colorbar
a=colorbar;
ylabel(a,'u')
colormap('jet')
title('u velocity of whole domain at Re=1, dt=0.001')
xlabel('x')
ylabel('y')

figure(4);
contourf(X,Y,p_final, 21, 'LineStyle', 'none')
colorbar
a=colorbar;
ylabel(a,'pressure')
colormap('jet')
title('pressure of whole domain at Re=1, dt=0.001')
xlabel('x')
ylabel('y')

figure(5);
hold on
quiver(X, Y, u_final, v_final, 5, 'k')
title('streamlines at Re=1, dt=0.001')
xlabel('x')
ylabel('y')

figure(6);
contourf(X,Y,cell_p{1,1}, 21, 'LineStyle', 'none')
colorbar
colormap('jet')
title('pressure of whole domain at Re=1, t=0.02')
xlabel('x')
ylabel('y')

figure(7);
contourf(X,Y,cell_p{1,2}, 21, 'LineStyle', 'none')
colorbar
a=colorbar;
ylabel(a,'pressure')
colormap('jet')
title('pressure of whole domain at Re=1, t=0.04')
xlabel('x')
ylabel('y')

figure(8);
contourf(X,Y,cell_p{1,3}, 21, 'LineStyle', 'none')
colorbar
a=colorbar;
ylabel(a,'pressure')
colormap('jet')
title('pressure of whole domain at Re=1, t=0.06')
xlabel('x')
ylabel('y')

figure(9);
contourf(X,Y,cell_p{1,4}, 21, 'LineStyle', 'none')
colorbar
a=colorbar;
ylabel(a,'pressure')
colormap('jet')
title('pressure of whole domain at Re=1, t=0.08')
xlabel('x')
ylabel('y')

figure(10);
contourf(X,Y,cell_p{1,5}, 21, 'LineStyle', 'none')
colorbar
a=colorbar;
ylabel(a,'pressure')
colormap('jet')
title('pressure of whole domain at Re=1, t=0.1')
xlabel('x')
ylabel('y')

figure(11);
hold on
quiver(X, Y, cell_u{1,1}, cell_v{1,1}, 5, 'r')
title('streamlines at Re=1, t=0.02')
xlabel('x')
ylabel('y')

figure(12);
hold on
quiver(X, Y, cell_u{1,2}, cell_v{1,2}, 5, 'r')
title('streamlines at Re=0.1, t=0.04')
xlabel('x')
ylabel('y')

figure(13);
hold on
quiver(X, Y, cell_u{1,3}, cell_v{1,3}, 5, 'r')
title('streamlines at Re=1, t=0.06')
xlabel('x')
ylabel('y')

figure(14);
hold on
quiver(X, Y, cell_u{1,4}, cell_v{1,4}, 5, 'r')
title('streamlines at Re=1, t=0.08')
xlabel('x')
ylabel('y')

figure(15);
hold on
quiver(X, Y, cell_u{1,5}, cell_v{1,5}, 5, 'r')
title('streamlines at Re=1, t=0.1')
xlabel('x')
ylabel('y')

figure(16);
hold on
plot(y2,u_05,'r--','LineWidth',2)
title('velocity profile at Re=1, t=0.1')
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k');

plot(y2,v_05,'k--','LineWidth',2)
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k');

Legend1={'u profile at x=0.5','v profile at x=0.5'};
legend(Legend1,'FontSize',15);
grid on
hold off
figure(17)
plot(y2,v_05,'k--','LineWidth',2)
title('v velocity profile at Re=1, t=0.1')
set(gca,'Fontsize',13,'linewidth',1.5)
xlabel('y','Fontsize',15,'FontWeight','bold','Color','k')
ylabel('velocity','Fontsize',15,'FontWeight','bold','Color','k');
Legend2={'v profile at x=0.5'};
legend(Legend2,'FontSize',15);
grid on