function cnf = riesz_sphere(cnf,N,dim,s,plotit)
%RIESZ_SPHERE
% cnf = riesz_sphere(cnf,N,dim,s,plotit)
% Returns a configuration obtained from applying the gradient descent to
% the given (or random) N-point collection on the unit sphere.
% Call without input arguments to use the defaults.
% 
% cnf -- pass your initial configuration as a matrix (dim)x(#of points); 
%   pass ZERO to draw from the Gaussian random distribution using your 
%   N,dim,s;
% N -- number of points in the random configuration to be generated
%   (ignored if an initial cnf is being passed);
% dim -- dimension of the ambient space; deduced from the first dimension
%   of the cnf matrix, if any;
% s -- exponent in the Riesz energy to be minimized; the default value is
%   5.0.
%   It is HIGHLY recommended to use either s=5.0 or s=0.5, as these are 
%   pre-coded, or to modify the source code. Otherwise you'll be using the 
%   Matlab's power function, which turns out to be not that great.
% plotit -- pass 'y' or 1, etc., to plot the produced configuration.
if ~exist('cnf','var')
    cnf = 1;
    N = 1000;
    dim = 3;
    s = 5.0;
    plotit = 1;
end
if cnf==0 || cnf==1
    fprintf( '\nStarting with a random point set.')
    cnf = randn(dim,N); 
else
    dim = size(cnf,1);
    N = size(cnf,2);
end
switch s
    case 5.0
        compute_weights = @(x) 1./x./x./x;
    case 2.5
        compute_weights = @(x) 1./x./sqrt(x)./sqrt(sqrt(x));
    case 0.5
        compute_weights = @(x) 1./sqrt(x)./sqrt(sqrt(x));
    otherwise
        compute_weights = @(x) sqrt(x).^(-s-1);     
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if s < dim
    k_value = N-1;  
    repel_steps = 100;
    repel_cycles = 5;
else
    k_value = 6 * dim;
    repel_steps = 200;
    repel_cycles = 10;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
fprintf( '\nWe will be minimizing the %3.2f-Riesz energy of %d points on the',s,N)
fprintf( '\n%d-dimensional unit sphere.\n\n', dim-1)
% double-check we're on the sphere:
cnf = cnf./sqrt(sum(cnf.*cnf,1));
% - this is faster than normc in Matlab 2016b

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
pbaspect([1 1 1])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',8)

for cycle=1:repel_cycles
    if s>=dim || cycle==1
        [IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
        IDX = IDX(:,2:end)';             % drop the trivial first column in IDX
        step = min(D(:,2));
    else
        if cycle<3
            [~, D] = knnsearch(cnf', cnf', 'k', 2);
            step = min(D(:,2));
        end
    end
    tic
    for iter=1:repel_steps
        cnf_repeated = reshape(repmat(cnf,k_value,1),dim,[]);
        knn_differences = cnf_repeated - cnf(:,IDX);
        knn_norms = sum(knn_differences.*knn_differences,1);
        riesz_weights = compute_weights(knn_norms);
        directions = bsxfun(@times,riesz_weights,knn_differences);
        directions = reshape(directions, dim, k_value, []);
        normals = reshape(sum(directions,2), dim, []);
%         normals = directions./sqrt(sum(directions.^2,1));
        tangents = normals-bsxfun(@times,sum(cnf.*normals,1),cnf); 
        % computing scalar products only works because we are on the unit
        % sphere
        tangents = tangents/max(sqrt(sum(tangents.^2,1)));
        cnf = cnf + tangents * step/3/cycle;%/iter
        cnf = cnf./sqrt(sum(cnf.*cnf,1));
    end
    toc
end
    
% [IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
% step = min(D(:,2));
msize = ceil(max(1, 22-3.5*log10(size(cnf,2)) ));
if dim==3 && (plotit=='y' || plotit=='Y' || plotit==1)
    colormap(winter)
    [x,y,z] = sphere(30);
    mesh(x,y,z)
    hold on
    plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',msize)
else
    if dim==2 && (plotit=='y' || plotit=='Y' || plotit==1)
        plot(cnf(1,:),cnf(2,:),'.k','MarkerSize',ceil(msize/2))
    end
end

% dlmwrite('cnf.out',cnf','delimiter','\t'); % ,'precision',3)