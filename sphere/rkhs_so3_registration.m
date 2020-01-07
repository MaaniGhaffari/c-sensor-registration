classdef rkhs_so3_registration < handle
    % RKHS_SO3_REGISTRATION class
    %   Author: Maani Ghaffari Jadidi and William Clark
    %   Date:   November 21, 2019
    % There is NO adaptive step size in this case!
    
    % PUBLIC PROPERTIES
    properties
        fixed;                  % fixed (target) point cloud
        moving;                 % moving (source) point cloud
        ell = 0.25;             % kernel characteristic length-scale
        sigma = 0.1;            % kernel signal variance (set as std)
        sp_threshold = 1e-3;    % kernel sparsification threshold
        c = 7;                  % so(3) inner product scale
        d = 7;                  % R^3 inner product scale
        color_scale = 1e-5;     % color space inner product scale
        MAX_ITER = 100;        % maximum number of iteration
        % The program stops if norm(omega)+norm(v) < eps
        eps = 5*1e-4;
        eps_2 = 1e-4;
        R = eye(3);             % initial orientation 
        T = zeros(3,1);         % initial translation
        omega;                  % so(3) part of twist
        v;                      % R^3 part of twist
        tform;                  % SE(3) tf
        iterations;             % number if iterations performed 
        ids;                    % nonzero indicies of coefficient matrix A
        min_step = 2*1e-1;      % minimum step size for integration
        step;                   % integration step
        length_x;               % target point cloud counts
        length_y;               % source point cloud counts
        cloud_x;                % target points
        cloud_y;                % source points
        CI;                     % matrix of color inner products
        A;                      % coefficient matrix
        K;
    end
    
    
    methods(Static)
        function CI = color_inner_product(Cx, Cz, scale)
            % Inner product matrix in color space.
            % CX and CZ are n-by-3 and m-by-3 color matrices of each point cloud,
            % respectively, where n and m are number of observations.
            % The resulting C matrix, is n-by-m and contains the inner product of all
            % combinations of points from X and Z.
            % scale is the parameter to scale the value of the inner product.
            
            
            %  Each column of C is n-by-1 vector of inner product of all Cx rows with a
            %  row of Cy. For example for the i-th column of C, we need
            % <Cx, cz_i'> = Cx * cz_i'. We can compute C all at once using Cx * Cz'.
            CI = scale * double(Cx) * double(Cz)';
        end
        
        function K = se_kernel(X, Z, l, s2)
            % Isotropic (same length-scale for all dimensions) squared-exponential kernel
            % X and Z are n-by-d and m-by-d inputs matrices, respectively, where n and
            % m are number of observations and d is the dimension of the input space.
            % The resulting kernel matrix, K, is n-by-m.
            % l and s2 are hyperparameters, i.e., length-scale and signal
            % variance, respectively.
            
            % number of observation for each input
            nx	= size(X,1);
            nz	= size(Z,1);
            
            % compute matrix of pair-wise squared distances
%             D2 = (sum((X.^2), 2) * ones(1,nz)) + (ones(nx, 1) * sum((Z.^2),2)') - ...
%                 2 * X * Z';
            
            DOT = X*Z'./((sum((X.^2), 2) * ones(1,nz)).*(ones(nx, 1) * sum((Z.^2),2)')); 
            D2 = real(acos(DOT));
            
            % Isotropic SE kernel
            K = s2 * exp(-D2.^2 / (2*l^2));
            % The new version
            K = K.*D2./sqrt(1-DOT.^2);
        end
        
        function Ti = tf_inv(R, t)
            % SE(3) inverse
            Ti = [R', -R' * t; 0, 0, 0, 1];
        end
        
        function d = dist_se3(R, t)
            d = norm(logm([R, t; 0, 0, 0, 1]),'fro');
        end
       
        function omega_hat = hat(omega)
            % The 'hat' function, R^3\to\Skew_3
            omega_hat = [0,         -omega(3),  omega(2);...
                        omega(3),   0,          -omega(1);...
                        -omega(2),  omega(1),   0];
        end
    end
    
    methods
        function obj = rkhs_se3_registration(varargin)
            % RKHS_SE3_REGISTRATION constructor
            if nargin == 2
                disp('Initial transformation is set.');
                obj.R = varargin{1};
                obj.T = varargin{2};
            elseif nargin > 0
                warning('The inputs are ignored!');
            end
        end
        
        function set_ptclouds(obj, target, source)
            if nargin == 3
                obj.fixed = target;
                obj.moving = source;
                obj.CI = obj.color_inner_product(target.Color, source.Color, obj.color_scale);
                obj.length_x = size(obj.fixed.Location,1);
                obj.length_y = size(obj.moving.Location,1);
                obj.cloud_x = double(obj.fixed.Location);
                obj.ell = 0.15;             % kernel characteristic length-scale
                obj.R = eye(3);
                obj.T = zeros(3,1);
            else
                error('Provide target and source point clouds as inputs');
            end
        end
        
        function compute_flow(obj)
            % Computes the Lie algebra transformation elements
            % twist = [omega; v]
            
            % Construct our kernel and coefficient
            obj.K = obj.se_kernel(obj.cloud_x, obj.cloud_y, obj.ell, obj.sigma^2);
            obj.K(obj.K < obj.sp_threshold) = 0;
            obj.A = sparse(double(obj.CI .* obj.K));
            
            % sum them up
            w = zeros(1,3); nu = w;
            used_ids = cell(obj.length_x,1);
            for i = 1:obj.length_x
                used_ids{i} = find(obj.A(i,:));
                Ai = nonzeros(obj.A(i,:))';
                cloud_yi = obj.cloud_y(used_ids{i},:);
                partial_omega = 1/obj.c * Ai * ...
                    cross(repmat(obj.cloud_x(i,:), length(used_ids{i}),1), cloud_yi ,2);
                partial_v = 1/obj.d * Ai * ...
                    (cloud_yi - obj.cloud_x(i,:));
                
                w = w + partial_omega;
                nu = nu + partial_v;
            end
            obj.ids = used_ids;
            obj.omega = w';
            obj.v = nu'; obj.v = 0*obj.v; % No translation
        end
        
        function compute_step_size(obj)
            % compute the integration step size
            obj.step = 0.1;
        end
        
        function align(obj)
            % Aligns two RGBD point clouds
            for k = 1:obj.MAX_ITER
                % construct omega and v
                % The point clouds are fixed (x) and moving (y)
                % current transformation
                obj.tform = affine3d(obj.tf_inv(obj.R,obj.T)');
                moved = pctransform(obj.moving, obj.tform);
                
                % extract point cloud information:
                obj.cloud_y = double(moved.Location);
                
                % compute omega and v
                obj.compute_flow;
                
                % compute step size for integrating the flow
                obj.compute_step_size;
                dt = obj.step;
                
                % Stop if the step size is too small
                if max(norm(obj.omega),norm(obj.v)) < obj.eps
                    break;
                end
                
                % Integrating
                th = norm(obj.omega); homega = obj.hat(obj.omega); 
                dR = eye(3) + (sin(dt * th) / th) * homega + ...
                    ((1 - cos(dt * th)) / th^2) * homega^2;
                dT = (dt * eye(3) + ...
                    (1-cos(dt * th)) / (th^2) * homega + ...
                    ((dt * th - sin(dt * th)) / th^3) * homega^2) * obj.v;
                R_new = obj.R * dR;
                T_new = obj.R * dT + obj.T;
                
                % Update the state
                obj.R = R_new;
                obj.T = T_new;
                
                % Our other break
                if obj.dist_se3(dR,dT) < obj.eps_2
                    break;
                end
                
                if k > 3
                    obj.ell = 0.15;
                end
                if k > 10
                    obj.ell = 0.10;
                end
                if k > 20
                    obj.ell = 0.05;
                end
                
                if k == 1
                   disp('   iter      grad_norm   eps      tf_dist     eps2') 
                end
                
                if mod(k,5) == 0
                    disp([k,max(norm(obj.omega),norm(obj.v)),obj.eps,obj.dist_se3(dR,dT),obj.eps_2]);
                end
            end
            
            obj.tform = affine3d(obj.tf_inv(obj.R, obj.T)');
            obj.iterations = k;
        end
    end
end