% Code to solve the 1D heat equation using Crank-Nicolson scheme
clc, clear,close

% Parameters
La=0.020875;    
dx=2; 
dt=0.1; 
nx=6;   
nt=11;

% Coefficients
b = La; 
c = b;        
a = 2*(1+La); 

% Initial and boundary conditions
Uo(1)=100;  Uo(2:nx-1)=0; Uo(nx)=50; Un(1)=100; Un(nx)=50;  

UUU(1,:)=Uo;

for k=2:nt    
    for ii=1:nx-2        
        if ii==1         
            d(ii)=c*Uo(ii)+2*(1-c)*Uo(ii+1)+b*Uo(ii+2)+c*Un(1);       
        elseif ii==nx-2           
            d(ii)=c*Uo(ii)+2*(1-c)*Uo(ii+1)+b*Uo(ii+2)+b*Un(nx);     
        else         
            d(ii)=c*Uo(ii)+2*(1-c)*Uo(ii+1)+b*Uo(ii+2);       
        end     
    end
    bb = b*ones(nx-3,1);
    cc = bb;
    aa = a*ones(nx-2,1);

    AA = diag(aa) + diag(-bb,1) + diag(-cc,-1);

    UU = AA\d';

    Un = [Un(1), UU', Un(nx)];

    UUU(k,:) = Un;

    Uo = Un;
end

UUU



