function cnf = riesz_torus(cnf,N,s,r,R,plotit,silent)
% RIESZ_TORUS
% cnf = riesz_torus(cnf,N,s,r,R)
% Returns a configuration obtained from performing the gradient descent on
% the given (or random) N-point collection on the torus with radii R and r,
% r <= R.
% 
% cnf -- pass your initial configuration as a matrix (dim)x(#of points)
%   here; 
%  a) call without input arguments to use all the default settings 
%   (cnf, 500, 6, 1, 3), cnf random as below;
%  b) pass ONE to draw from the uniform distribution IN ANGULAR measure 
%   (d\phi \times d\psi) using your values of N,s,r,R;
% N -- number of points in the random configuration to be generated
%   (ignored if an initial cnf is passed);
% s -- exponent in the Riesz energy to be minimized; default value is 6.0;
% r -- minor radius of the torus;
% R -- major radius of the torus.
% plotit -- pass 'y' or 1, etc., to plot the produced configuration.
% silent -- pass 'y' or 1, etc., to suppress output to console.
if ~exist('silent','var')
    silent = false;
end
if ~exist('plotit','var')
    plotit = 1;
end
if ~exist('s','var')
    s = 4.0;
end
if ~exist('r','var')
    r = 1.0;
end
if ~exist('R','var')
    R = 3.0;
end
if ~exist('cnf','var')
    cnf = 1;       
    N = 500;
end
k_value = 80;  
repel_steps = 100;
cycles = 8;
offset = 18;
torus = @(phi, theta,r,R) [ (R+r*cos(theta)).*cos(phi);...
                            (R+r*cos(theta)).*sin(phi);...
                            r*sin(theta)];
jtorus = @(phi, theta,r,R) [-(R+r*cos(theta)).*sin(phi); ...
                             (R+r*cos(theta)).*cos(phi); ...
                             zeros(size(theta));...          
                             -r*sin(theta).*cos(phi);...
                             -r*sin(theta).*sin(phi);...
                             r*cos(theta)];
if ~ismatrix(cnf) && (cnf==0 || cnf ==1)
    cnf = 2*pi*rand(2,N);
    cnf = torus(cnf(1,:),cnf(2,:),r,R);
else
    [~, N] = size(cnf);
end

if isscalar(cnf)
    if ~exist('silent','var') || ~silent
        fprintf( '\nStarting with a random point set.')
    end
    cnf = 2*pi*rand(2,N);
    cnf = torus(cnf(1,:),cnf(2,:),r,R);
else
    [dim, N] = size(cnf);
end

switch s
    case 4.0
        compute_weights = @(x) 1./x./x;
    case 2.0
        compute_weights = @(x) 1./x;
    case 0.5
        compute_weights = @(x) 1./sqrt(sqrt(x));
    otherwise
        compute_weights = @(x) sqrt(x).^(-s);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

if ~exist('silent','var') || ~silent
    fprintf( '\nWe will be minimizing the %f-Riesz energy of %d points on the',s,N)
    fprintf( '\n3-dimensional torus with radii R=%3.2f and r=%3.2f\n\n', R,r)
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
pbaspect([1 1 1])
colormap(spring)

syms phi theta;
x = (R+r*cos(theta))*cos(phi);
y = (R+r*cos(theta))*sin(phi);
z = r*sin(theta);
h=fsurf(x, y, z, [0 2*pi 0 2*pi], 'FaceAlpha',.9);
h.EdgeColor = 'none';
brighten(.9)
hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',8)
tic

for cycle=1:cycles
    if s>=dim || cycle<5
        [IDX, ~] = knnsearch(cnf', cnf', 'k', k_value+1);
        IDX = IDX(:,2:end)';             % drop the trivial first column in IDX
    end
    tic
    for iter=1:repel_steps*cycle
        cnf_repeated = reshape(repmat(cnf,k_value,1),dim,[]);
        knn_differences = cnf_repeated - cnf(:,IDX);
        knn_norms_squared = sum(knn_differences.*knn_differences,1);
        riesz_weights = s*compute_weights(knn_norms_squared)./knn_norms_squared;
        
        gradient = bsxfun(@times,riesz_weights,knn_differences);
        gradient = reshape(gradient, dim, k_value, []);
        directions = reshape(sum(gradient,2), dim, []);
        directions = directions/mean(sqrt(sum(directions.*directions,1)));
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %         
        [tangent_angles1, tangent_angles2] = ...
            torus_inversion(cnf(1,:),cnf(2,:),cnf(3,:),r,R);
        tangent_planes = jtorus(tangent_angles1,tangent_angles2,r,R);
        tangent_planes1 = tangent_planes(1:dim,:);
        tangent_planes2 = tangent_planes(dim+1:end,:);
%         the next four lines aren't this nice in higher dimensions
        tangent_planes1 = tangent_planes1./sqrt(sum(tangent_planes1.^2,1));
        tangent_planes2 = tangent_planes2./sqrt(sum(tangent_planes2.^2,1));
%         
        tangent_coefficients1 = sum(tangent_planes1 .* directions,1);
        tangent_coefficients2 = sum(tangent_planes2 .* directions,1);
        tangents =  bsxfun(@times,tangent_coefficients1,tangent_planes1)...
                  + bsxfun(@times,tangent_coefficients2,tangent_planes2);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %               
        step = sqrt(min(reshape(knn_norms_squared,k_value,[]),[],1));
%         tangents = tangents/max(sqrt(sum(tangents.*tangents,1)));
        cnf = cnf + tangents .* step/(offset+iter);
        cnfR =  [R*cnf(1:dim-1,:)./sqrt(sum(cnf(1:dim-1,:).^2,1)); zeros(1,N)];
        cnf = cnfR + r*(cnf-cnfR)./sqrt(sum((cnf-cnfR).^2,1));
    end
    if ~exist('silent','var') || ~silent
        toc
    end
end
msize = ceil(max(1, 22-3.5*log10(size(cnf,2)) ));
if dim==3 && exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
    plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',msize)
    pbaspect([1 1 1])
    daspect([1 1 1])
    set(gca, 'Clipping', 'off')
    axis vis3d
end
if ~usejava('desktop') && exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
    print(mfilename,'-dpdf','-r300','-bestfit')
end


% dlmwrite('cnf.out',cnf','delimiter','\t');



% % % % % % % % % % % % Sparse implementation
%         tic
%         spatangents = sparse(I(:), J(:), tangent_planes(:))\reshape(normals,[],1);
%         spatangents = bsxfun(@times, spatangents', reshape(tangent_planes,dim,[]));
%         spatangents = reshape(spatangents, dim, dim-1,[]);
%         spatangents = reshape(sum(spatangents, 2),dim,[]);
%         toc
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
