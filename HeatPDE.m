% Matlab Program 6: Heat Diffusion in one dimensional wire within the Fully
% Implicit Method
clear;
% Parameters to define the heat equation and the range in space and time
L = 1.; % Lenth of the wire
T =1.; % Final time
% Parameters needed to solve the equation within the fully implicit method
maxk = 1000; % Number of time steps
dt = 0.001;
n = 100; % Number of space steps
dx = 0.01;
cond = 4.; % Conductivity
b = cond*dt/(dx*dx); % Parameter of the method
% Initial temperature of the wire: a sinus.
for i = 1:(L/dx)+1
x(i) =(i-1)*dx;
u(i,1) =4-4*x(i);
end
% Temperature at the boundary (T=0)
for k=1:(T/dt)
time(k) = (k-1)*dt;
end
aa(1:n)=-b;
bb(1:n+1)=1.+2.*b;
cc(1:n)=-b;
M = diag(bb,0)+diag(aa,-1)+diag(cc,1);
M(1,n+1) = -b; M(n+1,1) = -b;
MM=inv(M);
% Implementation of the implicit method
for k=2:maxk % Time Loop
uu=u(1:n+1,k-1)+dt*(2*x(i)+1);
u(1:n+1,k)=MM*uu;
end
% Graphical representation of the temperature at different selected times
figure(1)
plot(x,u(:,1),'-',x,u(:,5),'-',x,u(:,10),'-',x,u(:,30),'-')
legend('t = 0.001','t = 0.005','t = 0.01','t = 0.03');
title('Temperature within the fully implicit method')
xlabel('X')
ylabel('T')
