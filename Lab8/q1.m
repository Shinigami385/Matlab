clear;
clc;

ic = @(x) sin(pi.*x);      %Initial Condition
bc_1 = @(t) 0;                                  %Boundary Condition 1
bc_2 = @(t) 0;                                  %Boundary Condition 2
a_x = 0;                                        %Range x 
b_x = 1;
a_t = 0;                                        %Range t
b_t = 0.1;
h = 0.1;
k = 0.005;                                 %k = h^2/4
c = 1;
exact_sol = @(x,t) exp(-(pi*pi*t)).*sin(pi*x);
%fprintf('\t\t FTCS\n');
[num_sol x t] = ftcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);

function [U x t] = ftcs(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2)
    
s = c*(k/(h^2));          %lambda = k/h^2;

x = [a_x:h:b_x];    %Discretize space
t = [a_t:k:b_t];    %Discretize time range
n = length(t);
m = length(x);
U = zeros(n,m);
U(1,:) = ic(x);
for i=2:n
    U(i,1) = bc_1(t(i));
    U(i,m) = bc_2(t(i));
    for j=2:m-1
        U(i,j) = s*U(i-1,j-1) + (1 - 2*s)*U(i-1,j) + s*U(i-1,j+1);
    end
end
end

function plot_end_time(exact_sol,num_sol,x,t)
figure;
[n m] = size(num_sol);
plot(x,num_sol(n,:),'r*');
hold on;
plot(x,exact_sol(x,t(n)),'b');
hold off;
legend('Numerical Solution','Exact Solution');
title('Plot of exact and Numerical Solution at Final time');


end

function surface_plot(exact_sol,num_sol,x,t)
    figure;
    subplot(1,2,1);
    mesh(x,t,num_sol);
    subplot(1,2,2);
    [X Y] = meshgrid(x,t);
    surf(X,Y,exact_sol(X,Y));
    title('Surface plot of Numerical Solution and Exact solution');
end
