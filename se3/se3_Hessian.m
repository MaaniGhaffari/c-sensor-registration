function H = se3_Hessian(fixed,moving,param)
% Function that computes the Hessian for the SE3 case
%
%   Author: William Clark
%   Date:   November 21, 2019

% fixed is a pointCloud
% moving is a pointCloud
% param = [ell, sigma]

ell = param(1); sigma = param(2);

% The kernel
k = @(x,y) sigma^2*exp(-norm(x-y)^2/(2*ell^2));

% Extract information from the clouds
length_x = size(fixed.Location,1);
length_y = size(moving.Location,1);
cloud_x = double(fixed.Location);
cloud_y = double(moving.Location);
color_x = double(fixed.Color); 
color_y = double(moving.Color);

H = zeros(6,6); % <- The Hessian
for i = 1:length_x
    for j = 1:length_y
        x = cloud_x(i,:); xC = color_x(i,:);
        z = cloud_y(j,:); zC = color_y(j,:);
        H = H + 1/ell^2*dot(xC,zC)*k(x,z)*[blockA(x,z,ell),blockC(x,z,ell)';blockC(x,z,ell),blockD(x,z,ell)];
    end
end

end


%% The blocks of the Hessian
function A = blockA(x,z,ell)
    y = cross(x,z);
    dot1 = x(2)*z(2)+x(3)*z(3); dot2 = x(1)*z(1)+x(3)*z(3);
    dot3 = x(1)*z(1)+x(2)*z(2);
    A11 = 1/ell^2*y(1)^2-dot1; A22 = 1/ell^2*y(2)^2-dot2; A33 = 1/ell^2*y(3)^2-dot3;
    A21 = 1/ell^2*y(1)*y(2)+1/2*(x(1)*z(2)+x(2)*z(1));
    A31 = 1/ell^2*y(1)*y(3)+1/2*(x(1)*z(3)+x(3)*z(1));
    A23 = 1/ell^2*y(2)*y(3)+1/2*(x(2)*z(3)+x(3)*z(2));
    A = [A11,A21,A31;A21,A22,A23;A31,A23,A33];
end

function C = blockC(x,z,ell)
    y = cross(x,z); v = z-x;
    C11 = 1/ell^2*y(1)*v(1); C22 = 1/ell^2*y(2)*v(2); C33 = 1/ell^2*y(3)*v(3);
    C21 = x(3)+1/ell^2*v(2)*y(1); C31 = -x(2)+1/ell^2*v(3)*y(1);
    C12 = -x(3)+1/ell^2*v(1)*y(2); C32 = x(1)+1/ell^2*v(3)*y(2);
    C13 = x(2)+1/ell^2*v(1)*y(3); C23 = -x(1)+1/ell^2*v(2)*y(3);
    C = [C11,C12,C13;C21,C22,C23;C31,C32,C33];
end

function D = blockD(x,z,ell)
    v = z-x;
    D11 = 1/ell^2*v(1)^2-1; D22 = 1/ell^2*v(2)^2-1; D33 = 1/ell^2*v(3)^2-1;
    D21 = 1/ell^2*v(1)*v(2); D31 = 1/ell^2*v(1)*v(3); D23 = 1/ell^2*v(2)*v(3);
    D = [D11,D21,D31;D21,D22,D23;D31,D23,D33];
end