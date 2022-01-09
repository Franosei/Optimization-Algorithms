%% Name: Francis Osei
%% Course: Optimization and Numerical Probability
%% Newtons Method with four variables
%% 20th December, 2021
clear;clc
syms X_1 X_2 X_3 X_4;
f = (X_1 + 10*X_2)^2 + 5*(X_3-X_4)^4 + (X_2-2*X_3)^4 + 10*(X_1-X_4)^4
e = 0.01;
%% Initial guess
xo = [3,-1,0,1] 
i = 1;
%% Computing the Gradient
dx_1 = diff(f, X_1);
dx_2 = diff(f, X_2);
dx_3 = diff(f, X_3);
dx_4 = diff(f, X_4);
g_func = [dx_1;dx_2;dx_3;dx_4];
g = [subs(dx_1,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)])...
subs(dx_2,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)])...
subs(dx_3,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)])...
subs(dx_4,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)])];
%% Computing the Hessian
dxx_1 =  diff(dx_1,X_1);
dxx_12 = diff(dx_1,X_2);
dxx_13 = diff(dx_1,X_3);
dxx_14 = diff(dx_1,X_4);
dxx_2 =  diff(dx_2,X_2);
dxx_23 = diff(dx_2,X_3);
dxx_24 = diff(dx_2,X_4);
dxx_3 =  diff(dx_3,X_3);
dxx_34 = diff(dx_3,X_4);
dxx_4 =  diff(dx_4,X_4);
dxx_111 = subs(dxx_1,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_121 = subs(dxx_12,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_131 = subs(dxx_13,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_141 = subs(dxx_14,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_222 = subs(dxx_2,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_232 = subs(dxx_23,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_242 = subs(dxx_24,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_333 = subs(dxx_3,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_343 = subs(dxx_34,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
dxx_444 = subs(dxx_4,[X_1,X_2,X_3,X_4],[xo(1),xo(2),xo(3),xo(4)]);
h = [dxx_111,dxx_121,dxx_131,dxx_141;dxx_121,dxx_222,dxx_232,dxx_242;...
[dxx_131,dxx_232,dxx_333,dxx_343];[dxx_141,dxx_242,dxx_343,dxx_444]];
h_inv = inv(h);
%% Objective function evaluated at the initial guess
f_x(i) = subs(f,[X_1,X_2,X_3,X_4], [xo(i),xo(i),xo(i),xo(i)]);

%% %%%%%%%%%%%%%%Implementing the newtions method%%%%%%%%%%%%%%%%%%%
while norm(g) > e
  w(1) = xo(1);x(1) = xo(2);y(1) = xo(3);z(1) = xo(4);
  I = [w(i),x(i),y(i),z(i)]'; %%%%Initial values
  w(i+1) = I(1) - h_inv(1,:)*g'; 
  x(i+1) = I(2) - h_inv(2,:)*g';
  y(i+1) = I(3) - h_inv(3,:)*g';
  z(i+1) = I(4) - h_inv(4,:)*g';
  i = i + 1; 
  %% updating the gradient
  g = [subs(dx_1,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)])...
  subs(dx_2,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)])...
  subs(dx_3,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)])...
  subs(dx_4,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)])];
  %% updating the hessian
  dxx_111 = subs(dxx_1,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]);
  dxx_121 = subs(dxx_12,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]); 
  dxx_131 = subs(dxx_13,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]);
  dxx_141 = subs(dxx_14,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]); 
  dxx_222 = subs(dxx_2,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]); 
  dxx_232 = subs(dxx_23,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]);
  dxx_242 = subs(dxx_24,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]);
  dxx_333 = subs(dxx_3,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]); 
  dxx_343 = subs(dxx_34,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]);
  dxx_444 = subs(dxx_4,[X_1,X_2,X_3,X_4],[w(i),x(i),y(i),z(i)]);
  h = [dxx_111,dxx_121,dxx_131,dxx_141;dxx_121,dxx_222,dxx_232,dxx_242;...
  [dxx_131,dxx_232,dxx_333,dxx_343];[dxx_141,dxx_242,dxx_343,dxx_444]];
  h_inv = inv(h);
  %% Evaluatinf the function for each point
  f_x(i) = subs(f,[X_1,X_2,X_3,X_4], [w(i),x(i),y(i),z(i)]); 
end
%% Representing the final outcome as a Table
itr = 1:i;
f_xk = f_x';
f_xk  = round(double(f_xk),10);
x_1 = w';x_2 = x';x_3 = y';x_4 = z';
iterations = itr';
T = table(x_1,x_2,x_3,x_4,f_xk,iterations)

