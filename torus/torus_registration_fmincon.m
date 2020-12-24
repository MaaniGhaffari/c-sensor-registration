% An example of continuous sensor registration on the torus. 
%
%   Author: Maani Ghaffari
%   Date:   December 24, 2020
clc; clear; close all

%% Parameters
numberX = 10;
numberZ = 15;
initial_disp = [1,1];
R = 3; r = 1; % <- Torus parameters

%% The clouds, it will be a 5 pointed star
% Define parameters
numberOfPoints = 5;
rotationAngle = 0;
xCenter = 0;
yCenter = 0;
% Determine the angles that the arm tips are at
theta = (0 : (numberOfPoints-1)/numberOfPoints*pi : (numberOfPoints-1)*pi) + rotationAngle;
% Define distance from the arm tip to the center of star.
amplitude = 1;
% Get x and y coordinates of the arm tips.
x = amplitude .* cos(theta) + xCenter;
y = amplitude .* sin(theta) + yCenter;

X = zeros(numberOfPoints * numberX,2); % The points
Z = zeros(numberOfPoints * 10,2);

% Go through the lines
for i = 1:numberOfPoints
    xx_part = linspace(x(i),x(i+1),numberX);
    xy_part = linspace(y(i),y(i+1),numberX);
    zx_part = linspace(x(i),x(i+1),numberZ) + rand(1,numberZ)/10 - 0.05;
    zy_part = linspace(y(i),y(i+1),numberZ) + rand(1,numberZ)/10 - 0.05;
    X((i-1)*numberX+1:i*numberX,:) = [xx_part',xy_part'];
    Z((i-1)*numberZ+1:i*numberZ,:) = [zx_part',zy_part'];
end

TX = makeTorus(X,R,r);

%% Perturb the problem
Z_wrong = Z + initial_disp;

%% Solve the registration
u0 = [0;0];
A = []; b = [];
Aeq = []; beq = [];
lb = [0;0]; ub = [2*pi;2*pi];
u = fmincon(@(x)objfunc_torus(x,X,Z_wrong),u0,A,b,Aeq,beq,lb,ub);

%% Set up the torus
th  = linspace(0,2*pi,64);
phi = linspace(0,2*pi,64);
[Phi,Th] = meshgrid(phi,th);
xT = (R+r.*cos(Th)).*cos(Phi);
yT = (R+r.*cos(Th)).*sin(Phi);
zT = r.*sin(Th);

%% The initial and final pictures
% Initial
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
TZ = makeTorus(Z_wrong,R,r);
surf(xT,yT,zT, 'EdgeColor', 'none'); alpha 0.35;
plot3(TX(:,1),TX(:,2),TX(:,3),'k.','LineWidth',2, 'MarkerSize', 14);
plot3(TZ(:,1),TZ(:,2),TZ(:,3),'r.','LineWidth',2, 'MarkerSize', 14);
axis equal tight; grid on; view(90,25);
hold off;
figuresize(18,18,'cm')
% title('Perturbed');
print(gcf,'intial_torus.png','-dpng','-r350');

% Final
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
TZ = makeTorus(Z_wrong - u',R,r);
surf(xT,yT,zT, 'EdgeColor', 'none'); alpha 0.35;
plot3(TX(:,1),TX(:,2),TX(:,3),'k.','LineWidth',2, 'MarkerSize', 14);
plot3(TZ(:,1),TZ(:,2),TZ(:,3),'r.','LineWidth',2, 'MarkerSize', 14);
axis equal tight; grid on; view(90,25);
hold off;
figuresize(18,18,'cm')
% title('Aligned');
print(gcf,'final_torus.png','-dpng','-r350');

%% Move to torus
function T = makeTorus(X,R,r)
    theta = X(:,1); phi = X(:,2);
    x = (R+r*cos(theta)).*cos(phi);
    y = (R+r*cos(theta)).*sin(phi);
    z = r*sin(theta);
    T = [x,y,z];
end

%% Objective function
function f = objfunc_torus(alpha,X,Z)
    f = 0; 
    XL = length(X); ZL = length(Z);
    Zt = Z - alpha';
    % Parameters
    ell = 0.25;
    sigma = 1;
    % Functions
    d2 = @(x,y) (cos(x)-cos(y))^2 + (sin(x)-sin(y))^2;
    dt = @(x,y) d2(x(1),y(1)) + d2(x(2),y(2)); % <- the distance on the torus
    k  = @(x,y) sigma^2*exp(-dt(x,y)/(2*ell^2));
    for i = 1:XL
        for j = 1:ZL
            f = f + k(X(i,:),Zt(j,:));
        end
    end
    f = -f; % maximize the objective
end