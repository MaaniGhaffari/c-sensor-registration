% A toy example of continuous sensor registration in R^2 (plane) with Lie
% group as SE(2).
%
%   Author: William Clark
%   Date:   November 21, 2019
clc; clear; close all

%% The fixed point cloud
x = linspace(-3,3,100); y = linspace(-3,3,100);
[X,Y] = meshgrid(x,y);
Z = peaks(100);
zmin = floor(min(Z(:)));
zmax = ceil(max(Z(:)));
zinc = (zmax-zmin)/40;
zlevs = zmin:zinc:zmax;

% The contour
figure; [C,~] = contour(X,Y,Z,zlevs);

% The color function
f = @(t) 1/16*[1,-1,0]*(t+7)+[0,1,0];

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
X = X'; X = [X,ones(length(X),1)]; XC = X;
for i = 1:length(XC)
    XC(i,:) = f(colX(i));
end
fixed = pointCloud(X,'Color',XC,'Intensity',colX');

%% The moving point cloud
x = linspace(-3,3,120); y = linspace(-3,3,120);
[X,Y] = meshgrid(x,y);
Z = peaks(120);
zmin = floor(min(Z(:)));
zmax = ceil(max(Z(:)));
zinc = (zmax-zmin)/43;
zlevs = zmin:zinc:zmax;

% The contour
figure(1); [C,h] = contour(X,Y,Z,zlevs);

% The color function
f = @(t) 1/16*[1,-1,0]*(t+7)+[0,1,0];

% Pull stuff out
u = length(C); loc = 1; Z = []; colZ = [];
while loc < u
    % The current level
    lvl = C(1,loc);
    % How far we need to go
    len = C(2,loc);
    % The next string
    Z = [Z,C(:,(loc+1):(loc+len))];
    colZ = [colZ,lvl*ones(1,len)];
    loc = loc + len + 1;
end
Z = Z'; Z = [Z,ones(length(Z),1)]; ZC = Z;
for i = 1:length(ZC)
    ZC(i,:) = f(colZ(i));
end
moving = pointCloud(Z,'Color',ZC,'Intensity',colZ');

%% Now we need to come up with an initial displacement
% Randomly perturb them
w = rand/2; v = rand(2,1);
A = [cos(w),-sin(w),0,v(1);sin(w),cos(w),0,v(2);0,0,1,0;0,0,0,1]';
moved = pctransform(moving, affine3d(A));

%% Now, recovering the pictures
rkhs_se2 = rkhs_se2_registration();
rkhs_se2.set_ptclouds(fixed,moved);
tic; rkhs_se2.align(); toc;
tform = rkhs_se2.tform;

%% Showing the results
aligned = pctransform(moved,tform);

%% The final errors
disp('Initial displacement:');
A = A'; T = tform.T'; B = [A(1:2,1:2),A(1:2,4);[0,0,1]];
disp(B);
disp('Algorithm finding:');
D = [T(1:2,1:2),T(1:2,4);[0,0,1]];
disp(D);
disp('Final displacement:');
disp(D*B);
disp('Frobenius norm of error:');
err = norm(logm(D*B),'fro');
disp(err);

%% Trying to plot something reasonable looking
% The X cloud
fixedXYZ = double(fixed.Location); fixedColor = double(fixed.Intensity);
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
scatter(fixedXYZ(:,1),fixedXYZ(:,2),1,fixedColor); 
axis equal tight; grid on;
% title('X');
figuresize(16,16,'cm')
print(gcf,'target_se2.png','-dpng','-r350');

% The Z cloud
movingXYZ = double(moving.Location); movingColor = double(moving.Intensity);
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
scatter(movingXYZ(:,1),movingXYZ(:,2),1,movingColor); 
axis equal tight; grid on;
% title('Z');
figuresize(16,16,'cm')
print(gcf,'source_se2.png','-dpng','-r350');

% The wrong pictures
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
% The X picture
scatter(fixedXYZ(:,1),fixedXYZ(:,2),1,fixedColor);
% The perturbed Z picture
movedXYZ = double(moved.Location); movedColor = double(moved.Intensity);
scatter(movedXYZ(:,1),movedXYZ(:,2),1,movedColor); 
axis equal tight; grid on;
% title('Perturbed');
figuresize(16,16,'cm')
print(gcf,'perturbed_se2.png','-dpng','-r350');

% The correct picture
figure; hold on; set(gca,'TickLabelInterpreter','latex', 'fontsize', 16); 
% The X picture
scatter(fixedXYZ(:,1),fixedXYZ(:,2),1,fixedColor);
% The moved Z picture
alignedXYZ = double(aligned.Location); alignedColor = double(aligned.Intensity);
scatter(alignedXYZ(:,1),alignedXYZ(:,2),1,alignedColor); 
axis equal tight; grid on;
% title('Aligned');
figuresize(16,16,'cm')
print(gcf,'aligned_se2.png','-dpng','-r350');