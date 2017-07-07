function cnf = riesz_shell(cnf,n,N,r,R,dim,s,plotit,analyzeit,silent)
%RIESZ_SHELL
% cnf = riesz_shell(cnf,n,N,r,R,dim,s,plotit,analyzeit,silent)
% Returns a configuration obtained from applying the gradient descent to
% the given (or random) 2*n+N-point collection inside a spherical shell.
% The inner and outer shell radii are defined by parameters r and R,
% respectively.
% Call without input arguments to use the defaults.
% 
% cnf -- pass your initial configuration as a matrix (dim)x(#of points); 
%   pass ZERO to draw from the random volume-uniform distribution using your 
%   values of n,N,r,R,dim,s;
% n -- number of points on the inner and outer boundaries (n points on
% each);
% N -- number of points in the interior of the shell;
%   (both n and N are ignored if an initial cnf is passed);
% r, R -- radii of the inner and outer boundary spheres; if a cnf array is
%   passed, but no radii specified, r and R are deduced as the min and max
%   vector norm of the columns of cnf;
% dim -- dimension of the ambient space; deduced from the first dimension
%   of the cnf matrix, if any;
% s -- exponent in the Riesz energy to be minimized; the default value is
%   5.0.
%   It is HIGHLY recommended to use one of the values {.5, 2.5, 5.0} for s,
%   as these are pre-coded, or to modify the source code. Otherwise you'll
%   be using the Matlab's power function, which turns out to be not that
%   great.
% plotit -- pass 'y' or 1, etc., to plot the produced configuration.
% analyzeit -- pass 'y' or 1, etc., to invoke pt_analyzer for the produced
% configuration.
% silent -- pass 'y' or 1, etc., to suppress output to console.
if ~exist('silent','var')
    silent = false;
end
if ~exist('analyzeit','var')
    analyzeit = 0;
end
if ~exist('plotit','var')
    plotit = 1;
end
if ~exist('s','var')
    s = 5.0;
end
if ~exist('dim','var')
    dim = 3;
end
if ~exist('cnf','var')
    cnf = 1;       
    n = 5000;
    N = 20000;
end
if isscalar(cnf)
    if ~exist('silent','var') || ~silent
        fprintf( '\nStarting with a random point set.')
    end
    if ~exist('r','var')
        r = 1;
    end
    if ~exist('R','var')
        R = 1.2;
    end
    cnf = randn(dim,N);
    % place cnf on the unit sphere
    cnf = cnf./sqrt(sum(cnf.*cnf,1));
    % - this is faster than normc in Matlab 2016b
    rads = (r^3 + (R^3-r^3)*rand(1,size(cnf,2))).^(1/3);
    % pushforward measure with the distribution function cubic in r and pdf
    % quadratic in r
    cnf = cnf .* rads;
    % cnf is now uniform in the volume measure
    fprintf( '\nGenerating surface distribution using the riesz_sphere routine.')
    rcnfsurf = r*riesz_sphere(0,n,dim,s,0,1);
    Rcnfsurf = R*riesz_sphere(0,n,dim,s,0,1);
    fprintf( '\nDone.')
    cnf = [cnf rcnfsurf Rcnfsurf];
else
    [dim, N] = size(cnf);
    if ~exist('r','var')
        r = min(sqrt(sum(cnf.*cnf,1)));
        n = 0;
    end
    if ~exist('R','var')
        R = max(sqrt(sum(cnf.*cnf,1)));
        n = 0;
    end
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
    repel_steps = 200;
    repel_cycles = 5;
else
    k_value = min(10 * dim, N-1);
    repel_steps = 50;
    repel_cycles = 8;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if ~exist('silent','var') || ~silent
    fprintf( '\nWe will be minimizing the %3.2f-Riesz energy of %d points in the',s,size(cnf,2))
    fprintf( '\n%d-dimensional shell with radii %3.2f and %3.2f.\n', dim-1, r, R)
    fprintf( 'This includes %d points on both inner and outer surfaces,\n', n)
    fprintf( 'and %d points in the interior.\n\n', N)
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',8)

for cycle=1:repel_cycles
    if s>=dim || cycle<4
        [IDX, ~] = knnsearch(cnf', cnf(:,1:N)', 'k', k_value+1);
        IDX = IDX(:,2:end)';             % drop the trivial first column in IDX
    end
    tic
    for iter=1:cycle*repel_steps
%       unitary = cnf./sqrt(sum(cnf.*cnf,1));
        cnf_repeated = reshape(repmat(cnf(:,1:N),k_value,1),dim,[]);
        knn_differences = cnf_repeated - cnf(:,IDX);
        knn_norms_squared = sum(knn_differences.*knn_differences,1);
        riesz_weights = compute_weights(knn_norms_squared);
        gradient = bsxfun(@times,riesz_weights,knn_differences);
        gradient = reshape(gradient, dim, k_value, []);
        directions = reshape(sum(gradient,2), dim, []);
        directions = directions/max(sqrt(sum(directions.*directions,1)));
        
        step = sqrt(min(reshape(knn_norms_squared,k_value,[]),[],1));
        
        cnf(:,1:N) = cnf(:,1:N) + directions .* step/iter;
        rads = sqrt(sum(cnf(:,1:N).*cnf(:,1:N),1));
        cnf(:,1:N) = cnf(:,1:N) ./ rads .* min(max(rads, r), R);
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
    daspect([1 1 1]);
    colormap(winter)
    [x,y,z] = sphere(30);
    X=r*x; Y=r*y; Z=r*z;
    mesh(X,Y,Z,'EdgeAlpha',.3,'FaceAlpha',1)
    hold on
    X=R*x; Y=R*y; Z=R*z;
    colormap(autumn)    
    mesh(X,Y,Z,'EdgeAlpha',.1,'FaceAlpha',.1)
    p1 = plot3(cnf(1,1:N),cnf(2,1:N),cnf(3,1:N),'.k','MarkerSize',msize);
    p2 = plot3(cnf(1,N+1:end),cnf(2,N+1:end),cnf(3,N+1:end),'ob','MarkerSize',.4*msize);
    leg=legend([p1; p2],'N interior nodes after redistribution (may be on the surface)',...
    '2*n initial surface nodes ');
    leg.Location = 'south';
    leg.FontSize = 12;
    axis vis3d
else
    if dim==2 && exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
        pbaspect([1 1 1]);
        plot(cnf(1,:),cnf(2,:),'.k','MarkerSize',ceil(msize/2))
    end
end
if ~usejava('desktop')
    print(mfilename,'-dpdf','-r300','-bestfit')
end
if dim==3 && exist('analyzeit','var') && (analyzeit=='y' || analyzeit=='Y' || analyzeit==1)
    in_d = @(x,y,z) min(max(sqrt(x.*x + y.*y + z.*z), r), R) == sqrt(x.*x + y.*y + z.*z);
    pt_analyzer(cnf, in_d)
end

% dlmwrite('cnf.out',cnf','delimiter','\t','precision',10);