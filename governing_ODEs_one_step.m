clear all; close all;

% System parameters used in paper 
p = [0.5, 0.0035, 0.019, 9.81, 0.25, 0.016];  % [A1, k, a1, g, A2, a2]
global A1 A2 k a1 a2 g
A1 = p(1);
k = p(2);
a1 = p(3);
g = p(4);
A2 = p(5);
a2 = p(6);

% load twotankdata
% U=u(1:end)
% U_val=u(2001:end);
% Y_actual=y(1:end);
% Y_val=y(2001:end);

U=zeros(3000,1);
U(1:1000)=1;
U(1001:2000)=2;
U(2001:3000)=3;
%Y=zeros(2000,1);

dt = 0.2;                 % Time step
num_steps = length(U);    % Total number of time steps from U

% x0 = [0.1730; 0.2439];  % Initial steady state values found for upper and lower tank water levels
x0 = [0.0017; 0.0024];
myoptions = odeset('NonNegative',[1:2]);
Y(1) = 0;
Ty = 0;
for i = 1:num_steps
    u = U(i);  % Using input from loaded data after initial steady state found
    % u=10;     % Used to find initial steady state
    [tout,xout] = ode45(@(t,x) odefun(t,x,u),[(i-1)*dt i*dt],x0,myoptions);
    Y = [Y;xout(end,2)];
    Ty = [Ty;tout(end)];
    x0 = xout(end,:)';
end

time = (0:num_steps-1)' * dt;  % Time vector

figure;
subplot(2, 1, 1);
plot(time, U, 'r-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Input Pump Voltage U (V)');
title('Step Input for Pump Voltage');
grid on;

subplot(2, 1, 2);
plot(time, Y(2:end), 'b-', 'LineWidth', 1.5);
xlabel('Time [s]');
ylabel('Water Level in Lower Tank Y (m)');
title('Response of Lower Tank Water Level to Step Input');
grid on;

Y(1) = Y(2);

U_val = U(1501:end);
Y_val = Y(1502:end);
U = U(1:1500);
Y = Y(2:1501);
save('OneStepTwoTanksMatlab.mat','U','Y','U_val','Y_val')

function xdot = odefun(t,x,u)
global A1 A2 k a1 a2 g
    dx1 = (1/A1) * (k * u - a1 * sqrt(2 * g * max(0,x(1)))); % d/dt x1
    dx2 = (1/A2) * (a1 * sqrt(2 * g * x(1)) - a2 * sqrt(2 * g * max(0,x(2)))); % d/dt x2
    xdot = [dx1;dx2];
end