function cnf = riesz_shell(cnf,varargin)
%RIESZ_SHELL
% cnf = riesz_shell(cnf,'NAME1',VALUE1,...,'NAMEN',VALUEN)
% Returns a configuration obtained from applying the gradient descent to
% the given (or random) 2*n+N-point collection inside a spherical shell.
% The inner and outer shell radii are defined by parameters r and R,
% respectively. The 2n points on the boundaries are held fixed.
% Call without input arguments to use the defaults.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%
% cnf -- pass your initial configuration as a matrix (dim)x(#of points);
%   pass ZERO to draw from the random vol-uniform distribution using your
%   M, N, r, R, dim, s;
% Optional argument name/value pairs:
% Name          Value
%
% 'M'           number of points held on the inner and outer boundaries (M
%               points on each). If cnf is passed, deduced as the mean of
%               the number of points with the min and max norms; default: 0
% 'N'           number of points to be generated in the interior of the
%               shell (ignored if a cnf is passed); default: 20,000
% 'r', 'R'      radii of the inner and outer boundary spheres; if a cnf
%               array is passed, r and R are deduced as the min and max
%               vector norm of the columns of cnf; defaults: 1.0 and 1.2
% 'dim'         dimension of the ambient space; deduced from the first
%               dimension of the cnf matrix, if given; default: 3
% 'plotit'      pass 'y' or 1, true, etc., to plot the produced
%               configuration; default: true
% 'silent'      pass 'y' or 1, true, etc., to suppress output to console;
%               default: false
% 'analyzeit'	pass 'y' or 1, true, etc., to invoke f_analyzer on the
%               produced configuration; default: false
% 's'           the exponent used in the Riesz kernel;
%               default: 4.0
%   It is HIGHLY recommended to use s from {0.5, 2.0, 4.0}, as these are
%   pre-coded, or to modify the source code. Otherwise you'll be using the
%   Matlab's power function, which turns out to be not that great.


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Initialize variables
pnames = {'M'   'N'     'r'  'R' 'dim' 'plotit' 'silent' 'analyzeit' 's' };
dflts =  {0  20000   1.0  1.2   3    true    false    false       4.0 };
[M, N, r, R, dim, plotit, silent, analyzeit, s, ~] =...
    internal.stats.parseArgs(pnames, dflts, varargin{:});
if ~exist('cnf','var') || isscalar(cnf)
    if ~silent
        fprintf( '\nStarting with a random point set.')
    end
    cnf = randn(dim,N);
    % place cnf on the unit sphere:
    cnf = cnf./sqrt(sum(cnf.*cnf,1));
    % - this is faster than normc in Matlab 2016b
    rads = (r^3 + (R^3-r^3)*rand(1,size(cnf,2))).^(1/3);
    % pushforward measure with the distribution function cubic in r and pdf
    % quadratic in r
    cnf = cnf .* rads;
    % cnf is now uniform in the volume measure
    if M > 0
        fprintf( '\nGenerating surface distribution using the riesz_sphere routine.')
        rcnfsurf = r*riesz_sphere(0,'N', M,'dim', dim,'s', s, 'plotit', 0, 'silent',1);
        Rcnfsurf = R*riesz_sphere(0,'N', M,'dim', dim,'s', s, 'plotit', 0, 'silent',1);
        fprintf( '\nDone.')
        cnf = [cnf rcnfsurf Rcnfsurf];
    end
else
    [ dim, N ]= size(cnf);
    cnfrs = sqrt(sum(cnf.*cnf,1));
    r = min(cnfrs);
    R = max(cnfrs);
    wtbig = 0.999;
    wtsmall = 1 - wtbig;
    num1 = cnfrs < (r * wtbig + R* wtsmall);
    num2 = cnfrs > (r * wtsmall + R* wtbig);
    M = floor(mean([sum(num1) sum(num2)]));
    N = N - 2 * M;
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
%% Flow cycle parameters
offset = 18;
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
%% Talk
if ~silent
    fprintf( '\nWe will be minimizing the %3.2f-Riesz energy of %d points in the',s,size(cnf,2))
    fprintf( '\n%d-dimensional shell with radii %3.2f and %3.2f.\n', dim-1, r, R)
    fprintf( 'This includes %d points on both inner and outer surfaces,\n', M)
    fprintf( 'and %d points in the interior.\n\n', N)
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
% plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',8)
%% Main

for cycle=1:repel_cycles
    if s>=dim || cycle<4
        [IDX, ~] = knnsearch(cnf', cnf(:,1:N)', 'k', k_value+1);
        IDX = IDX(:,2:end)';             % drop the trivial first column in IDX
    end
    tic
    for iter=1:cycle*repel_steps
        cnf_repeated = reshape(repmat(cnf(:,1:N),k_value,1),dim,[]);
        knn_differences = cnf_repeated - cnf(:,IDX);
        knn_norms_squared = sum(knn_differences.*knn_differences,1);
        riesz_weights = compute_riesz(knn_norms_squared)./knn_norms_squared;
        gradient = bsxfun(@times,riesz_weights,knn_differences);
        gradient = reshape(gradient, dim, k_value, []);
        directions = reshape(sum(gradient,2), dim, []);
        directions = directions./sqrt(sum(directions.*directions,1));
        
        step = sqrt(min(reshape(knn_norms_squared,k_value,[]),[],1));
        
        cnf(:,1:N) = cnf(:,1:N) + directions .* step/(offset + iter);
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
    colormap(winter)
    [x,y,z] = sphere(30);
    X=r*x; Y=r*y; Z=r*z;
    mesh(X,Y,Z,'EdgeAlpha',.3,'FaceAlpha',1)
    hold on
    X=R*x; Y=R*y; Z=R*z;
    colormap(autumn)
    mesh(X,Y,Z,'EdgeAlpha',.1,'FaceAlpha',.1)
    p1 = plot3(cnf(1,1:N),cnf(2,1:N),cnf(3,1:N),'.k','MarkerSize',msize);
    s1 = N + " interior nodes after redistribution (may be on the surface)";
    s2 = "2*" + M + " fixed surface nodes";
    if M > 0
        p2 = plot3(cnf(1,N+1:end),cnf(2,N+1:end),cnf(3,N+1:end),'ob','MarkerSize',.4*msize);
        leg=legend([p1; p2],s1,s2);
    else
        leg=legend(p1,s1);
    end
    leg.Location = 'south';
    leg.FontSize = 12;
    pbaspect([1 1 1])
    daspect([1 1 1])
    set(gca, 'Clipping', 'off')
    axis vis3d
else
    if dim==2 && exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
        pbaspect([1 1 1]);
        plot(cnf(1,:),cnf(2,:),'.k','MarkerSize',ceil(msize/2))
    end
end
if ~usejava('desktop') && exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
    print(mfilename,'-dpdf','-r300','-bestfit')
end
if dim==3 && exist('analyzeit','var') && (analyzeit=='y' || analyzeit=='Y' || analyzeit==1)
    in_d = @(x,y,z) min(max(sqrt(x.*x + y.*y + z.*z), r), R) == sqrt(x.*x + y.*y + z.*z);
    f_analyzer(cnf, in_d)
end

% dlmwrite('cnf.out',cnf','delimiter','\t','precision',10);
