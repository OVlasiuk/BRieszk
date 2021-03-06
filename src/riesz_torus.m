function out = riesz_torus(cnf, varargin)
% RIESZ_TORUS
% cnf = riesz_torus(cnf,'NAME1',VALUE1,...,'NAMEN',VALUEN)
% Returns a configuration obtained from performing the gradient descent on
% the given (or random) N-point collection on the torus with radii R and r,
% r <= R.
% Call without input arguments to use the defaults.
% Large outputs without output arguments (assignment) are suppressed.
% 
% cnf -- pass your initial configuration as a matrix (dim)x(#of points); 
%   pass ZERO to draw a random distribution using your N, dim, r, R, s;
% Optional argument name/value pairs:
% Name          Value
%
% 'N'           number of points in the random configuration to be generated
%               (ignored if an initial cnf is passed); default: 500
% 'dim'         dimension of the ambient space; deduced from the first dimension
%               of the cnf matrix, if given; default: 3
% 'r'           minor radius of the torus; default: 1.0
% 'R'           major radius of the torus; default: 3.0
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
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% Initialize variables 
pnames={'N' , 'dim' , 'r' , 'R' , 'plotit' , 'silent' , 'saveit' , 'saveplot' , 's'};
dflts={1000 , 3     , 1.0 , 3.0 , true     , true     , false    , false      , 4.0};
[N , dim , r , R , plotit , silent , saveit , saveplot , s, ~] =...
     internal.stats.parseArgs(pnames, dflts, varargin{:});
%% Maps
torus = @(phi, theta,r,R) [ (R+r*cos(theta)).*cos(phi);...
                            (R+r*cos(theta)).*sin(phi);...
                            r*sin(theta)];
jtorus = @(phi, theta,r,R) [-(R+r*cos(theta)).*sin(phi); ...
                             (R+r*cos(theta)).*cos(phi); ...
                             zeros(size(theta));...          
                             -r*sin(theta).*cos(phi);...
                             -r*sin(theta).*sin(phi);...
                             r*cos(theta)];
 
%% Configuration 
if ~exist('cnf','var') || isscalar(cnf)
    if ~silent
        fprintf( '\nStarting with a random point set.')
    end
	cnf = 2*pi*rand(2,N);
    cnf = torus(cnf(1,:),cnf(2,:),r,R);   
else
    [ dim, N ]= size(cnf);
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
%% Flow cycle parameters
offset = 18;
if s < dim
    k_value = N-1;  
    repel_steps = 100;
    cycles = 8;
else
    k_value = min(6 * dim, N-1);
    repel_steps = 50;
    cycles = 6;
end

%% Talk  % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

if ~exist('silent','var') || ~silent
    fprintf( '\nWe will be minimizing the %f-Riesz energy of %d points on the',s,N)
    fprintf( '\n3-dimensional torus with radii R=%3.2f and r=%3.2f\n\n', R,r)
end

%% Main % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
pbaspect([1 1 1])

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
            f_torinv(cnf(1,:),cnf(2,:),cnf(3,:),r,R);
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
if dim==3 && exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
    pplot(cnf)
end

%% Output and saving 
if nargout()>0 || N<50
    out=cnf;
else
    out = "Refusing to output into console; see helpstring.";
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
