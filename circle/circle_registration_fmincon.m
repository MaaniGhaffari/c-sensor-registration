% An example of continuous sensor registration on a circle.
%
%   Author: Maani Ghaffari
%   Date:   December 24, 2020

clc; clear; close all

%% The clouds
X = linspace(0,pi/2,10);
Z = linspace(pi/2,pi,12);

u0 = 0;
A = []; b = [];
Aeq = []; beq = [];
lb = 0; ub = 2*pi;
u = fmincon(@(x)objfunc_circle(x,X,Z),u0,A,b,Aeq,beq,lb,ub);

%% Save the answers
% Initial
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',2);
plot(cos(X),sin(X),'k*',cos(Z),sin(Z),'r*','LineWidth',2);
grid on, axis equal tight
hold off;
figuresize(9,9,'cm')
print(gcf,'initial_circle.png','-dpng','-r350');

% Final
Zf = Z - u;
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',2);
plot(cos(X),sin(X),'k*',cos(Zf),sin(Zf),'r*','LineWidth',2);
grid on, axis equal tight
hold off;
figuresize(9,9,'cm')
print(gcf,'final_circle.png','-dpng','-r350');

%% Objective function
function f = objfunc_circle(alpha,X,Z)
    f = 0; 
    XL = length(X); ZL = length(Z);
    Zt = Z - alpha;
    % Parameters
    ell = 0.25;
    sigma = 1;
    % Functions
    d2 = @(x,y) (cos(x)-cos(y))^2 + (sin(x)-sin(y))^2;
    k = @(x,y) sigma^2*exp(-d2(x,y)/(2*ell^2));
    for i = 1:XL
        for j = 1:ZL
            f = f + k(X(i),Zt(j));
        end
    end
    f = -f; % maximize the objective
end