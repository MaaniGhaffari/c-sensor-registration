% An example of continuous sensor registration on a circle. We are simply
% integrating the gradient in this case. 
%
%   Author: William Clark
%   Date:   November 21, 2019

clc; clear; close all

%% Parameters
ell = 0.25;
sigma = 1;

%% Functions
d2 = @(x,y) (cos(x)-cos(y))^2 + (sin(x)-sin(y))^2;
k = @(x,y) sigma^2*exp(-d2(x,y)/(2*ell^2));
grad_part = @(x,y) 1/ell^2*k(x,y)*sin(y-x);

%% The clouds
X = linspace(0,pi/2,10);
Z = linspace(pi/2,pi,12);

%% Perform the registration problem
diff_ode = @(t,u)odegrad(X,Z,u,grad_part);
[T,U] = ode45(diff_ode,[0,0.1],0);

%% Animate everything
v = VideoWriter('circle.avi','Uncompressed AVI');
open(v);
figure('Renderer', 'painters', 'Position', [0 0 720 480])
for i = 1:length(T)
    u_current = U(i);
    Z_current = Z - u_current;
    clf;
    hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
    plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',2);
    plot(cos(X),sin(X),'k*',cos(Z_current),sin(Z_current),'r*','LineWidth',2);
    grid on, axis equal tight
    set(gca,'visible','off')
    frame = getframe(gcf, [0 0 720 480]);
    writeVideo(v,frame);
    pause(1e-6);
end
close(v);
close all

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
Zf = Z - U(end);
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
plot(cos(linspace(0,2*pi)),sin(linspace(0,2*pi)),'k','LineWidth',2);
plot(cos(X),sin(X),'k*',cos(Zf),sin(Zf),'r*','LineWidth',2);
grid on, axis equal tight
hold off;
figuresize(9,9,'cm')
print(gcf,'final_circle.png','-dpng','-r350');


%% The differential equation
function df = odegrad(X,Z,alpha,grad_part)
    u = 0; XL = length(X); ZL = length(Z);
    Zt = Z - alpha;
    for i = 1:XL
        for j = 1:ZL
            u = u + grad_part(X(i),Zt(j));
        end
    end
    df = u;
end