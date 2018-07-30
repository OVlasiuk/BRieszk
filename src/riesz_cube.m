function out = riesz_cube(cnf,varargin)
%RIESZ_SPHERE
% cnf = riesz_sphere(cnf,'NAME1',VALUE1,...,'NAMEN',VALUEN)
% Returns a configuration obtained from applying the gradient descent to
% the given (or random) N-point collection on the unit cube.
% Call without input arguments to use the defaults.
% Large outputs without output arguments (assignment) are suppressed.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% 
% cnf -- pass your initial configuration as a matrix (dim)x(#of points); 
%   pass ZERO to draw from the uniform random distribution using your 
%   N, dim, s;
% Optional argument name/value pairs:
% Name          Value
%
% 'N'           number of points in the random configuration to be generated
%               (ignored if an initial cnf is passed); default: 1000
% 'dim'         dimension of the ambient space; deduced from the first dimension
%               of the cnf matrix, if given; default: 3
% 'plotit'      pass 'y' or 1, true, etc., to plot the produced 
%               configuration; default: true
% 'silent'      pass 'y' or 1, true, etc., to suppress output to console;
%               default: true
% 'saveplot'    pass 'y' or 1, true, etc., to save the plot;
%               default: false
% 'saveit'      pass 'string' to save output into a file named 'string';
%               default: false
% 's'           the exponent used in the Riesz kernel;
%               default: 4.0
%   It is HIGHLY recommended to use s from {0.5, 2.0, 4.0}, as these are 
%   pre-coded, or to modify the source code. Otherwise you'll be using the 
%   Matlab's power function, which turns out to be not that great.

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Initialize variables
pnames={'N' , 'dim' , 'plotit' , 'silent' , 'saveit' , 'saveplot' , 's'};
dflts={1000 , 3     , true     , true     , false    , false      , 4.0};
[N, dim, plotit, silent, saveit, saveplot, s, ~] =...
     internal.stats.parseArgs(pnames, dflts, varargin{:});
if ~exist('cnf','var') || isscalar(cnf)
    if ~silent
        fprintf( '\nStarting with a random point set.')
    end
    cnf = rand(dim,N); 
else
    [dim, N] = size(cnf);
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
offset = 10;
if s < dim
    k_value = N-1;  
    repel_steps = 100;
    repel_cycles = 5;
else
    k_value = min(10 * dim, N-1);
    repel_steps = 20;
    repel_cycles = 50;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Talk
if ~exist('silent','var') || ~silent
    fprintf( '\nWe will be minimizing the %3.2f-Riesz energy of %d points on the',s,N)
    fprintf( '\n%d-dimensional unit cube.\n\n', dim-1)
end

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%% Main

for cycle=1:repel_cycles
    if s>=dim || mod(cycle,3)==1
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
        %         tangents = gradient-bsxfun(@times,sum(cnf.*gradient,1),cnf); 
        % computing scalar products like so only works because we are on 
        % the unit sphere
        gradient = gradient./(sqrt(sum(gradient.^2,1)));
        step = sqrt(min(reshape(knn_norms_squared,k_value,[]),[],1));
        %         alpha = min(1, step./sqrt(sum(tangents.*tangents,1)));
        cnf = cnf + gradient .* step/(offset+iter);
        cnf(cnf>1) = 2-cnf(cnf>1);
        cnf(cnf<0) = -cnf(cnf<0);
    end
    if ~exist('silent','var') || ~silent
        toc
    end
end

msize = ceil(max(1, 22-3.5*log10(size(cnf,2)) ));
if dim==3 && exist('plotit','var') && (plotit)
    pbaspect([1 1 1]);
    colormap(winter)
    plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',msize)
    axis vis3d
else
    if dim==2 && exist('plotit','var') &&(plotit=='y' || plotit=='Y' || plotit==1)
        pbaspect([1 1 1]);
        plot(cnf(1,:),cnf(2,:),'.k','MarkerSize',ceil(msize/2))
    end
end
if ~usejava('desktop') && (plotit=='y' || plotit=='Y' || plotit==1 || plotit == true)
    print(mfilename,'-dpdf','-r300','-bestfit')
end

%% Output and saving 
if nargout()>0 || N<50
    out=cnf;
else
    out = "Refusing to output into console; see docstring by calling 'help "...
        +mfilename+"'";
end

if saveplot
    savehere = string(mfilename) + "_s_" + string(s) + "_N_" + string(N);
    disp("Saving figure to file: " + savehere)
    print(char(savehere),'-dpdf','-r300','-bestfit')
end


if isstr(saveit)
    disp("Saving configuration to file: " + saveit)
    dlmwrite(saveit,cnf','delimiter','\t','precision',10)
end
