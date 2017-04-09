function cnf = riesz_torus(cnf,N,s,r,R)
% RIESZ_TORUS
% cnf = riesz_torus(cnf,N,s,r,R)
% Returns a configuration obtained from performing the gradient descent on
% the given (or random) N-point collection on the torus with radii R and r,
% r <= R.
% 
% cnf -- pass your initial configuration as a matrix (dim)x(#of points)
%   here; 
%  a) pass ZERO to use all default settings (cnf, 10000, 6, 1, 3), cnf as 
%   below;
%  b) pass ONE to draw from the uniform distribution IN ANGULAR measure 
%   (d\phi \times d\psi) using your values of N,s,r,R;
% N -- number of points in the random configuration to be generated
%   (ignored if an initial cnf is passed);
% s -- exponent in the Riesz energy to be minimized; default value is 6.0;
% r -- minor radius of the torus;
% R -- major radius of the torus.
dim = 3;
k_value = 20;  
repel_steps = 100;
F = 3*repel_steps^0.6;
cycles = 40;
torus = @(phi, theta,r,R) [ (R+r*cos(theta)).*cos(phi);...
                            (R+r*cos(theta)).*sin(phi);...
                            r*sin(theta)];
jtorus = @(phi, theta,r,R) [-(R+r*cos(theta)).*sin(phi); ...
                             (R+r*cos(theta)).*cos(phi); ...
                             zeros(size(theta));...          
                             -r*sin(theta).*cos(phi);...
                             -r*sin(theta).*sin(phi);...
                             r*cos(theta)];
if cnf==0 || cnf==1
    s = cnf*s + (1-cnf)*6.0;
    N = cnf*N + (1-cnf)*10000;
    r = cnf*r + (1-cnf)*1.0;
    R = cnf*R + (1-cnf)*3.0;
    cnf = 2*pi*rand(2,N);
    cnf = torus(cnf(1,:),cnf(2,:),r,R);
else
    [~, N] = size(cnf);
end

fprintf( '\nWe will be minimizing the %f-Riesz energy of %d points on the',s,N)
fprintf( '\n3-dimensional torus with radii R=%3.2f and r=%3.2f\n\n', R,r)

% row indices in the sparse matrix with Jacobians:
I= repmat(reshape(1:N*dim,dim,[]), dim-1, 1);
% column indices in the sparse matrix with Jacobians:
J = repmat(1:(dim-1)*N, dim, 1);
tangent_coefficients = zeros(dim-1, N);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
pbaspect([1 1 1])
colormap(spring)

syms phi theta;
x = (R+r*cos(theta))*cos(phi);
y = (R+r*cos(theta))*sin(phi);
z = r*sin(theta);
h=fsurf(x, y, z, [0 2*pi 0 2*pi]);
h.EdgeColor = 'none';
brighten(.9)
hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',8)
tic
cnf1=cnf;

for cycle=1:cycles
    [IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
    IDX = IDX(:,2:end)';             % drop the trivial first column in IDX
    step = min(D(:,2));
%      tic
    fprintf( '\nNow in cycle %d.',cycle)
    for iter=1:repel_steps
        cnf_repeated = reshape(repmat(cnf,k_value,1),dim,[]);
        knn_differences = cnf_repeated - cnf(:,IDX);
        knn_norms = sqrt(sum(knn_differences.^2,1));
        riesz_weights = knn_norms.^(-s-1);
        directions = bsxfun(@times,riesz_weights,knn_differences);
        directions = reshape(directions, dim, k_value, []);
        normals = reshape(sum(directions,2), dim, []);
        [tangent_angles1, tangent_angles2] = ...
            torus_inversion(cnf(1,:),cnf(2,:),cnf(3,:),r,R);
        tangent_planes = jtorus(tangent_angles1,tangent_angles2,r,R);
        tangent_planes1 = tangent_planes(1:dim,:);
        tangent_planes2 = tangent_planes(dim+1:end,:);
%         the next four lines aren't this nice in higher dimensions
        tangent_planes1 = tangent_planes1./sqrt(sum(tangent_planes1.^2,1));
        tangent_planes2 = tangent_planes2./sqrt(sum(tangent_planes2.^2,1));
%         
        tangent_coefficients1 = sum(tangent_planes1 .* normals,1);
        tangent_coefficients2 = sum(tangent_planes2 .* normals,1);
        tangents =  bsxfun(@times,tangent_coefficients1,tangent_planes1)...
                  + bsxfun(@times,tangent_coefficients2,tangent_planes2);
%               
        tangents = tangents./sqrt(sum(tangents.^2,1));
        cnf = cnf + tangents * step/F;
        cnfR =  [R*cnf(1:dim-1,:)./sqrt(sum(cnf(1:dim-1,:).^2,1)); zeros(1,N)];
        cnf = cnfR + r*(cnf-cnfR)./sqrt(sum((cnf-cnfR).^2,1));
    end
%     toc 
end
mean_shift = mean(sqrt(sum((cnf-cnf1).^2,1)))
toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
% [IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
% step = min(D(:,2));
plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',8)
% plot3( cnf1(1,:), cnf1(2,:), cnf1(3,:),'.r','MarkerSize',7)
% dlmwrite('cnf.out',cnf','delimiter','\t');



% % % % % % % % % % % % Sparse implementation
%         tic
%         spatangents = sparse(I(:), J(:), tangent_planes(:))\reshape(normals,[],1);
%         spatangents = bsxfun(@times, spatangents', reshape(tangent_planes,dim,[]));
%         spatangents = reshape(spatangents, dim, dim-1,[]);
%         spatangents = reshape(sum(spatangents, 2),dim,[]);
%         toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 