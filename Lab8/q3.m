clear;
clc;

ic = @(x) 2*min(x,1-x);      %Initial Condition
bc_1 = @(t) 0;                                  %Boundary Condition 1
bc_2 = @(t) 0;                                  %Boundary Condition 2
a_x = 0;                                        %Range x 
b_x = 1;
a_t = 0;                                        %Range t
b_t = 0.1;
h = 0.1;
k = 0.01;                                 %k = h^2/4
c = 1;
exact_sol = @(x,t) exp(-(pi*pi*t)).*sin(pi*x);

fprintf('\t\tCrank Nicolson\n');
[num_sol x t] = Crank_Nicolson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2);
plot_end_time(exact_sol,num_sol,x,t);
surface_plot(exact_sol,num_sol,x,t);

function [U x t] = Crank_Nicolson(c,a_x,b_x,a_t,b_t,h,k,ic,bc_1,bc_2)
s = c*(k/(h^2));            %lambda = k/h^2;

x = [a_x:h:b_x];            %Discretize space
t = [a_t:k:b_t];            %Discretize time range
n = length(t);
m = length(x);
U = zeros(n,m);
U(1,:) = ic(x);
%U(t,x)
%(2I + sB)U(j,:) = (2I-sB)U(j-1,:)
B = zeros(m,m);
B(1,1) = 2;
B(1,2) = -1;
B(m,m-1) = -1;
B(m,m) = 2;
for i = 2:m-1
    B(i,i-1) = -1;
    B(i,i) = 2;
    B(i,i+1) = -1;
end
U = U';
for j = 2:n
   U(1,j-1) = bc_1(t(j-1));
   U(m,j-1) = bc_2(t(j-1));
   U(:,j) = (2*eye(m,m) + s*B)\(2*eye(m,m) - s*B)*U(:,j-1); 
end
U(1,n) = bc_1(t(n));
U(m,n) = bc_2(t(n));
U = U';
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

