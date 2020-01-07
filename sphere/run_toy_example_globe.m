% A toy example of continuous sensor registration on the sphere with Lie
% group as SO(3).
%
%   Author: William Clark
%   Date:   November 21, 2019
clc; clear; close all

%% Import the globe
load topo;
x = linspace(0,2*pi,360); y = linspace(0,pi,180);
[X,Y] = meshgrid(x,y);
% plot the globe
figure; 
[C,~] = contour(X,Y,topo);

% The color function
f = @(t) 1/14000*[1,-1,0]*(t+8000)+[0,1,0];

% Pull stuff out
u = length(C); loc = 1; X = []; colX = [];
while loc < u
    % The current level
    lvl = C(1,loc);
    % How far we need to go
    len = C(2,loc);
    % The next string
    X = [X,C(:,(loc+1):(loc+len))];
    colX = [colX,lvl*ones(1,len)];
    loc = loc + len + 1;
end
X = X';
X = [X,ones(length(X),1)]; XC = X;
for i = 1:length(XC)
    XC(i,:) = f(colX(i));
end
% X is currently on the plane, we need to transform it to the sphere
f_sphere = @(x,y)[sin(pi-y).*cos(x),sin(pi-y).*sin(x),cos(pi-y)];
X_sphere = f_sphere(X(:,1),X(:,2));
fixed = pointCloud(X_sphere,'Color',XC,'Intensity',colX');
fixed_down = pcdownsample(fixed,'gridAverage',0.05);

%% The moving point cloud
moving = fixed; % <- Let's make them the same
moving_down = pcdownsample(moving,'gridAverage',0.05);

%% Now we need to come up with an initial displacement
% Randomly perturb them
w = rand(3,1)/3; W = [0,-w(3),w(2);w(3),0,-w(1);-w(2),w(1),0];
A = expm(W); A_tform = [A,[0;0;0];0,0,0,1]';
moved = pctransform(moving, affine3d(A_tform));
moved_down = pctransform(moving_down, affine3d(A_tform));

%% Now, recovering the pictures
rkhs_s2 = rkhs_so3_registration();
rkhs_s2.set_ptclouds(fixed_down,moved_down);
tic; rkhs_s2.align(); toc;
tform = rkhs_s2.tform;

%% Showing the results
aligned = pctransform(moved,tform);

%% The final errors
disp('Initial displacement:');
T = tform.T';
disp(A);
disp('Algorithm finding:');
D = T(1:3,1:3);
disp(D);
disp('Final displacement:');
disp(D*A);
disp('Frobenius norm of error:');
err = norm(logm(D*A),'fro');
disp(err);

%% Trying to plot something reasonable looking
% The X cloud
fixedXYZ = double(fixed.Location); fixedColor = double(fixed.Intensity);
figure; scatter3(fixedXYZ(:,1),fixedXYZ(:,2),fixedXYZ(:,3),1,fixedColor);
% title('Reference'); 
axis tight square off;
figuresize(16,16,'cm')
print(gcf,'target.png','-dpng','-r350');

% The wrong pictures
figure;
% The X picture
scatter3(fixedXYZ(:,1),fixedXYZ(:,2),fixedXYZ(:,3),1,fixedColor);
hold on;
% The perturbed Z picture
movedXYZ = double(moved.Location); movedColor = double(moved.Intensity);
scatter3(movedXYZ(:,1),movedXYZ(:,2),movedXYZ(:,3),1,movedColor);
% title('Perturbed'); 
axis tight square off;
figuresize(16,16,'cm')
print(gcf,'perturbed.png','-dpng','-r350');

% The correct picture
figure;
% The X picture
scatter3(fixedXYZ(:,1),fixedXYZ(:,2),fixedXYZ(:,3),1,fixedColor);
hold on;
% The moved Z picture
alignedXYZ = double(aligned.Location); alignedColor = double(aligned.Intensity);
scatter3(alignedXYZ(:,1),alignedXYZ(:,2),alignedXYZ(:,3),1,alignedColor);
% title('Aligned'); 
axis tight square off;
figuresize(16,16,'cm')
print(gcf,'aligned.png','-dpng','-r350');