function [guess]=Iteration_initial_condition_Blasius_equation

%parameters
unk=0.0001;%delta_eta
delta_eta=unk;
length_eta=8;
n_eta=length_eta/delta_eta;

%iteration for an initial condition for y''
guess=1;
for i=1:1000
    [y_1,y_2,y_3]=Runge_kutta_method(guess,unk);
    ksi=1e-10;
    err=abs(y_2(n_eta+1)-1);
    if err>ksi
        guess=guess-0.3*err;
    end
end
end