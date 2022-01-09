%% Name: Francis Osei
%% Course: Optimization and Numerical Probability
%% Projected Gradient Method with four variables
%% Date: 20th December,2021
%%
format short;
clc;clear;
syms x_1 x_2 x_3 x_4;
%% Ojective function 
f = (x_1+x_2)^2+(x_3+x_4)^4+3*(x_1-2)^2+ (2*x_3+1)^2+2*(x_4-0.5)^2+1
p_f = 0.5*[(x_1-x_2),(x_2-x_1),(x_3-x_4),(x_4-x_3)]
%% initial guess
xo = [1 -1 1 -1]
e = 0.01;
%% 
i = 1;
%% gradient
g_func = gradient(f);
grad = (subs(g_func,[x_1,x_2,x_3,x_4],[xo(1),xo(2),xo(3),xo(4)]));
xk = xo
%x_k = subs(p_f,[x_1,x_2,x_3,x_4],[xo(1),xo(2),xo(3),xo(4)]);
%% projected gradient descent
while norm(grad)>e
    w(i) = xk(1);x(i) = xk(2);y(i) = xk(3);z(i) = xk(4);
    I = [w(i),x(i),y(i),z(i)];
    x1 = I(1) - 0.2*grad(1); 
    x2 = I(2) - 0.2*grad(2);
    x3 = I(3) - 0.2*grad(3);
    x4 = I(4) - 0.2*grad(4);
    J = [x1, x2, x3, x4];
    %% Updating the gradient
    xk = subs(p_f,[x_1,x_2,x_3,x_4],[J(1),J(2),J(3),J(4)]);
    grad = (subs(g_func,[x_1,x_2,x_3,x_4],[xk(1),xk(2),xk(3),xk(4)]));
    f_x(i) = subs(f,[x_1,x_2,x_3,x_4], [xk(1),xk(2),xk(3),xk(4)]);
    i = i + 1;
end
%% Representing the final outcome as a Table
itr = 1:i-1;
f_xk = f_x';
f_xk  = round(double(f_xk),10);
x_1 = w';x_2 = x';x_3 = y';x_4 = z';
iterations = itr';
T = table(x_1,x_2,x_3,x_4,f_xk,iterations)