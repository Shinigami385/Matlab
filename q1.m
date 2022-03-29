
clear;
clc;
f = @(l1,l2,l3,y) 0.5*[ l2+l3 l3-l1 l2-l1 ; l3-l2 l1+l3 l1-l2 ; l2-l3 l1-l3 l1+l2]*y ;
sol = @(l1,l2,l3,x) [exp(l2*x) + exp(l3*x) - exp(l1*x); exp(l1*x) + exp(l3*x) - exp(l2*x);exp(l2*x) + exp(l1*x) - exp(l3*x)];

l1 = -1;
l2 = 0;
l3 = 1;
fprintf("Lamda1 = %d , Lamda2 = %d , Lamda3 = %d \n",l1,l2,l3);
%fprintf("%2s %15s %15s %15s %15s \n",'N','EM','RK4','Error(EM)','Error(RK4)');
fprintf(" N \t\t EM \t\t\t \t RK4 \t\t\t Error(EM) \t\t\t Error(RK4)\n");
for N = [10 20 40 80]
    % EM
    y=[1;1;1];
    for i = 1:N
        y = y + (1/N)*f(l1,l2,l3,y);
    end

    % RK4
    rk4=[1;1;1];
    for i = 1:N
        k1 = (1/N)*f(l1,l2,l3,rk4);
        k2 = (1/N)*f(l1,l2,l3,rk4+(k1/2));
        k3 = (1/N)*f(l1,l2,l3,rk4+(k2/2));
        k4 = (1/N)*f(l1,l2,l3,rk4+k3);
        rk4 = rk4 + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
%fprintf("%2.1f [%2.6f;%2.6f;%2.6f] %15f %15f %15f \n",N,y,rk4,abs(y-sol(l1,l2,l3,1)),abs(rk4-sol(l1,l2,l3,1)));
fprintf("%2.1f [%2.6f;%2.6f;%2.6f] [%2.6f;%2.6f;%2.6f] [%2.6f;%2.6f;%2.6f] [%2.6f;%2.6f;%2.6f] \n", ...
    N,y,rk4,abs(y-sol(l1,l2,l3,1)),abs(rk4-sol(l1,l2,l3,1)));
end

l1 = 0;
l2 = -1;
l3 = -10;
fprintf("Lamda1 = %d , Lamda2 = %d , Lamda3 = %d \n",l1,l2,l3);
%fprintf("%2s %15s %15s %15s %15s \n",'N','EM','RK4','Error(EM)','Error(RK4)');
fprintf(" N \t\t EM \t\t\t \t RK4 \t\t\t Error(EM) \t\t\t Error(RK4)\n");
for N = [10 20 40 80]
    % EM
    y=[1;1;1];
    for i = 1:N
        y = y + (1/N)*f(l1,l2,l3,y);
    end

    % RK4
    rk4=[1;1;1];
    for i = 1:N
        k1 = (1/N)*f(l1,l2,l3,rk4);
        k2 = (1/N)*f(l1,l2,l3,rk4+(k1/2));
        k3 = (1/N)*f(l1,l2,l3,rk4+(k2/2));
        k4 = (1/N)*f(l1,l2,l3,rk4+k3);
        rk4 = rk4 + (k1 + 2*k2 + 2*k3 + k4)/6;
    end
    fprintf("%2f \n",N);
    %fprintf(y,rk4,abs(y-sol(l1,l2,l3,1)),abs(y-sol(l1,l2,l3,1)));
    fprintf("%2.1f [%2.6f;%2.6f;%2.6f] [%2.6f;%2.6f;%2.6f] [%2.6f;%2.6f;%2.6f] [%2.6f;%2.6f;%2.6f] \n", ...
        N,y,rk4,abs(y-sol(l1,l2,l3,1)),abs(rk4-sol(l1,l2,l3,1)));
end

