function [y_1,y_2,y_3]=Runge_kutta_method(guess,unk)

%parameters
delta_eta=unk;
length_eta=8;
n_eta=length_eta/delta_eta;

%construct answer matrix
y_1=zeros(1,n_eta+1);%y
y_2=zeros(1,n_eta+1);%first order differential of y
y_3=zeros(1,n_eta+1);%second order differential of y

%initial value
y_1(1)=0;
y_2(1)=0;
y_3(1)=guess;

%runge kutta method
for i=1:n_eta
    k_11=y_2(i);%y_1
    k_12=y_2(i)+0.5*delta_eta*k_11;
    k_13=y_2(i)+0.5*delta_eta*k_12;
    k_14=y_2(i)+delta_eta*k_13;
    k_21=y_3(i);%y_2
    k_22=y_3(i)+0.5*delta_eta*k_21;
    k_23=y_3(i)+0.5*delta_eta*k_22;
    k_24=y_3(i)+delta_eta*k_23;
    k_31=-0.5*y_1(i)*y_3(i);%y_3
    k_32=-0.5*(y_1(i)+0.5*delta_eta*k_31)*(y_3(i)+0.5*delta_eta*k_31);
    k_33=-0.5*(y_1(i)+0.5*delta_eta*k_32)*(y_3(i)+0.5*delta_eta*k_32);
    k_34=-0.5*(y_1(i)+0.5*delta_eta*k_33)*(y_3(i)+0.5*delta_eta*k_33);
    y_1(i+1)=y_1(i)+(delta_eta/6)*(k_11+2*k_12+2*k_13+k_14);
    y_2(i+1)=y_2(i)+(delta_eta/6)*(k_21+2*k_22+2*k_23+k_24);
    y_3(i+1)=y_3(i)+(delta_eta/6)*(k_31+2*k_32+2*k_33+k_34);
end
end
