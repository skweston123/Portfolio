% Devon Doud, Nathan Sanders, Sam Weston
% MSU Spring, 2019
% EMEC 303 Project 1

clear; clc
close all
 for k = [60,73.8,90]
% Set step size and length of time in seconds
h  = 1;
Lt = 6000;
N  = round(Lt/h);

% Preallocation 
height      = zeros(1,N);    
den         = zeros(1,N);
t           = zeros(1,N);
phi         = zeros(1,N);
force_r     = zeros(1,N);
force_theta = zeros(1,N);
m           = zeros(1,N);

% Define initial conditions and constants
t_m    = 1.1045;        % Thrust multiplier
t(1)   = 0;             % Initial time
den(1) = rho(0);        % Starting density
r(1)   = 6.378e6;       % Radius of earth
vr(1)  = 0;             % Inital radial velocity
th(1)  = 0;             % Initial angle theta
w(1)   = 7.29e-5;       % Angular velocity of earth
phi(1) = deg2rad(k); % Launch angle relative to horizontal
height(1) = 0;          % Initial heigh above ground

% Call force function to determine additional initial conditions
[force_r(1),force_theta(1),m(1)] =...
    force(t(1),den(1),phi(1),vr(1),w(1),r(1),t_m);

% Store dependent variables into y
y = [r,vr,th,w];

% ODEs describing motion of rocket
f = @(t,y,force_r,force_theta,m)...
    [y(2), force_r/m + y(1)*y(4)^2,... 
    y(4),(force_theta/m - 2*y(2)*y(4))/y(1)];

for i = 1:N-1
    
    % Update height and density 
    height(i+1) = y(i,1)-r(1);
    den(i+1) = rho(height(i));
    
    % Update time
    t(i+1)   = t(i) + h;
    
    % Update phi based on height
    if height(i) <= 275000
        phi(i+1) = sqrt(phi(1)-phi(1)*(height(i)/275000)^2);
    else
        phi(i+1) = 0;
    end
    
    % Call force function to update force in r and theta directions
    [force_r(i+1),force_theta(i+1),m(i+1)] = ...
        force(t(i),den(i),phi(i),y(i,2),y(i,4),y(i,1),t_m);
    
    % Euler update for dependent variables
    y(i+1,:) = y(i,:) + h*f(t(i),y(i,:),force_r(i),force_theta(i),m(i));

end
if k == 60
    y1=y;
elseif k == 73.8
    y2=y;
else
    y3=y;
end

end
% Plot results
% figure
% polarplot((y(:,3)),y(:,1))
% hold on
% polarplot(0:.01:2*pi,ones(size(0:.01:2*pi))*r(1))
% title('Rocket Trajectory in Polar Coordinates')
% word1={'Radius (m)'};
% text(0,3000000,word1)
% 
% figure
% subplot(2,2,1)
% plot(height,phi)
% title ('Rocket Angle Relative to Horizontal')
% ylabel('Angle (rad)')
% xlabel('Height (m)')
% 
% subplot(2,2,2)
% plot(t,height)
% title ('Rocket Height vs Time')
% ylabel('Height (m)')
% xlabel('Time (s)')
% 
% subplot(2,2,3)
% plot(t,y(:,2))
% title ('Radial Velocity vs Time')
% ylabel('Speed (m/s)')
% xlabel('Time (s)')
% 
% subplot(2,2,4)
% plot(t,y(:,4).*y(:,1))
% title ('Velocity in the Theta Direction vs Time')
% ylabel('Speed (m/s)')
% xlabel('Time (s)')

polarplot((y1(:,3)),y1(:,1))
hold on
polarplot((y2(:,3)),y2(:,1))
polarplot((y3(:,3)),y3(:,1))
title('Rocket Trajectory in Polar Coordinates')
word1={'Radius (m)'};
text(0,3000000,word1)
polarplot(0:.01:2*pi,ones(size(0:.01:2*pi))*r(1))
legend('60 degrees','73.8 degrees','90 degrees')