function cnf = riesz_sphere(cnf,N,dim,s,plotit,silent)
%RIESZ_SPHERE
% cnf = riesz_sphere(cnf,N,dim,s,plotit,silent)
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
% s -- the exponent used in the Riesz kernel;
%   It is HIGHLY recommended to use either s=4.0 or s=0.5, as these are 
%   pre-coded, or to modify the source code. Otherwise you'll be using the 
%   Matlab's power function, which turns out to be not that great.
% plotit -- pass 'y' or 1, etc., to plot the produced configuration.
% silent -- pass 'y' or 1, etc., to suppress output to console.
if ~exist('cnf','var')
    cnf = 1;
    N = 1000;
    dim = 3;
    s = 5.0;
    plotit = 1;
    silent = false;
end
if isscalar(cnf)
    if ~exist('silent','var') || ~silent
        fprintf( '\nStarting with a random point set.')
    end
    if ~exist('N','var')
        N = 1000;
    end
    if ~exist('dim','var')
        dim = 3;
    end
    cnf = randn(dim,N); 
else
    [dim, N] = size(cnf);
end
if ~exist('s', 'var')
    s = 4.0;
end
switch s
    case 4.0
        compute_riesz = @(x) 1./x./x;
    case 2.0
        compute_riesz = @(x) 1./x;
    case 0.5
        compute_riesz = @(x) 1./sqrt(sqrt(x));
    otherwise
        compute_riesz = @(x) sqrt(x).^(-s);
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if s < dim
    k_value = N-1;  
    repel_steps = 100;
    repel_cycles = 5;
else
    k_value = min(6 * dim, N-1);
    repel_steps = 50;
    repel_cycles = 10;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if ~exist('silent','var') || ~silent
    fprintf( '\nWe will be minimizing the %3.2f-Riesz energy of %d points on the',s,N)
    fprintf( '\n%d-dimensional unit sphere.\n\n', dim-1)
end
% double-check we're on the sphere:
cnf = cnf./sqrt(sum(cnf.*cnf,1));
% - this is faster than normc in Matlab 2016b

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',8)

for cycle=1:repel_cycles
    if s>=dim || cycle==1
        [IDX, ~] = knnsearch(cnf', cnf', 'k', k_value+1);
        IDX = IDX(:,2:end)';             % drop the trivial first column in IDX
    end
    tic
    for iter=1:repel_steps*cycle
        cnf_repeated = reshape(repmat(cnf,k_value,1),dim,[]);
        knn_differences = cnf_repeated - cnf(:,IDX);
        knn_norms_squared = sum(knn_differences.*knn_differences,1);
        riesz_weights = compute_riesz(knn_norms_squared)./knn_norms_squared;
        gradient = bsxfun(@times,riesz_weights,knn_differences);
        gradient = reshape(gradient, dim, k_value, []);
        gradient = reshape(sum(gradient,2), dim, []);
        tangents = gradient-bsxfun(@times,sum(cnf.*gradient,1),cnf); 
        % computing scalar products like so only works because we are on 
        % the unit sphere
        tangents = tangents/max(sqrt(sum(tangents.^2,1)));
        step = sqrt(min(reshape(knn_norms_squared,k_value,[]),[],1));
        cnf = cnf + tangents .* step/(3+iter);
        cnf = cnf./sqrt(sum(cnf.*cnf,1));
    end
    if ~exist('silent','var') || ~silent
        toc
    end
end

% [IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
% step = min(D(:,2));
msize = ceil(max(1, 22-3.5*log10(size(cnf,2)) ));
if dim==3 && exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
    pbaspect([1 1 1]);
    colormap(winter)
    [x,y,z] = sphere(30);
    mesh(x,y,z,'EdgeAlpha',.3)
    hold on
    plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',msize)
    axis vis3d
else
    if dim==2 && exist('plotit','var') &&(plotit=='y' || plotit=='Y' || plotit==1)
        pbaspect([1 1 1]);
        plot(cnf(1,:),cnf(2,:),'.k','MarkerSize',ceil(msize/2))
    end
end

fprintf('Compute the full Riesz energy of this pointset? [y/N]\n')
inp = input('','s');
if ~isempty(inp) && ((inp=='y') || (inp=='Y') || (inp=='1'))
    en =    (cnf(1,:)-cnf(1,:)').*(cnf(1,:)-cnf(1,:)') +...
                (cnf(2,:)-cnf(2,:)').*(cnf(2,:)-cnf(2,:)') +...
                (cnf(3,:)-cnf(3,:)').*(cnf(3,:)-cnf(3,:)');
    en = compute_riesz(en);
    en = sum(sum(en(isfinite(en))));
    fprintf('The %3.2f-Riesz energy is \t %10.6f\n', s,en)
end


% dlmwrite('cnf.out',cnf','delimiter','\t'); % ,'precision',3)