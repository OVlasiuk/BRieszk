% function [vorFig, triFig] = riesz_surf(cnf, surfF, gradF, )
%RIESZ_SURF
%  
% INPUT:
% cnf -- 3x(num_pts), the node set to be processed
% gradF -- a function handle to evaluate the gradient to the surface at
%   a point (x,y,z);
% OUTPUT:
% vorFig -- handle to the figure with the surface Voronoi diagram;
% triFig -- handle to the surface triangulation with vertices at cnf.
% Both are returned with the .Visible attribute set to 'off'.
%

% surfF = @(x) x(1).^2 .*(x(1).^2 - 5) + x(2).^2 .*(x(2).^2 - 5) +...
%                 x(3).^2 .*(x(3).^2 - 5) + 11;

% pnames = { 'jitter' 'pullback' 'A'     's'     'histogram' 'bins' 'offset' 'instats'};
% dflts =  { 0            []      100.0   4.0     false       200    18       true};
% [jitter, pullbackF, A, s, htrue, bins, offset, instats, ~] =...
%      internal.stats.parseArgs(pnames, dflts, varargin{:});


% EXAMPLE of a density defined by (the absolute value of) the Gaussian
% curvature:
% Goldman, R. (2005). Curvature formulas for implicit curves and surfaces.
% Computer Aided Geometric Design, 22(7), 632–658. 
% doi:10.1016/j.cagd.2005.06.005

% clear surfF;
cnf = dlmread('../cnf40k_2.txt'3')';

surfF = @(x) x(1,:).^2 .*(x(1,:).^2 - 5) + x(2,:).^2 .*(x(2,:).^2 - 5) +...
    x(3,:).^2 .*(x(3,:).^2 - 5) + 11;
gradF = @(x) [...
    2*x(1,:).*(x(1,:).*x(1,:) - 5) + 2*x(1,:).^3;
    2*x(2,:).*(x(2,:).*x(2,:) - 5) + 2*x(2,:).^3;
    2*x(3,:).*(x(3,:).*x(3,:) - 5) + 2*x(3,:).^3];

complementedlaplacianF = @(x) [...
    (12*x(2,:).*x(2,:) - 10) .* (12*x(3,:).*x(3,:) - 10);
    (12*x(1,:).*x(1,:) - 10) .* (12*x(3,:).*x(3,:) - 10);
    (12*x(1,:).*x(1,:) - 10) .* (12*x(2,:).*x(2,:) - 10)...
    ];

squaredgradientF = @(x) gradF(x).^2;

densityF = @(x)  sum(complementedlaplacianF(x).* squaredgradientF(x), 1)./ ...
    sum(squaredgradientF(x), 1).^2;

gg = @(x,y,z) [2*x*(x^2 - 5) + 2*x^3
2*y*(y^2 - 5) + 2*y^3
2*z*(z^2 - 5) + 2*z^3];

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

N = 40000;
% cnf = cnfinit(N, surfF);
% cnf = dlmread('../sconf40k.txt')';
N_moving = size(cnf,2);
% N = size(cnf,2);

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% plot3(x(1,:),x(2,:),x(3,:),'.b')
% pbaspect([1 1 1])
% daspect([1 1 1])
% axis vis3d
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

s = 4.0;
repel_steps = 1000;
dim = 3;
offset = 100;
k_value = 30;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
switch s
    case 4.0
        compute_riesz = @(x) 1./x./x;
        compute_weights = @(x) x.*x.*x.*x;
    case 2.0
        compute_riesz = @(x) 1./x;
        compute_weights = @(x) x.*x;
    case 0.5
        compute_riesz = @(x) 1./sqrt(sqrt(x));
        compute_weights = @(x) sqrt(x);
    otherwise
        compute_riesz = @(x) sqrt(x).^(-s);
        compute_weights = @(x) x.^s;
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
tic
ngrad = gradF(cnf);
%% Main loop
for iter=1:repel_steps
    if mod(iter,10) == 1
        [IDX, ~] = knnsearch(cnf', cnf(:,1:N_moving)', 'k', k_value+1);
        IDX = IDX(:,2:end)';
    end
%% Vectors from nearest neighbors    
    cnf_repeated = reshape(repmat(cnf(:,1:N_moving),k_value,1),dim,[]);
    knn_cnf = cnf(:,IDX);
    knn_differences = cnf_repeated - knn_cnf;    
    knn_norms_squared = sum(knn_differences.*knn_differences,1); 
%% Weights using radial density
    riesz_weights = compute_riesz(knn_norms_squared);
    if isa(densityF,'function_handle')
        knn_density =  abs(densityF(knn_cnf)) + 1;   
        weights = s* riesz_weights ./ knn_norms_squared ./ knn_density;
    else
        weights = s*riesz_weights./knn_norms_squared;
    end
%% Sum up over the nearest neighbors    
    gradient = bsxfun(@times,weights,knn_differences);
    gradient = reshape(gradient, dim, k_value, []);
    gradient = reshape(sum(gradient,2), dim, []);
%     
    surfnormals = ngrad./sqrt(sum(ngrad.*ngrad,1));
    tangentgrad = gradient - surfnormals .* sum(gradient.*surfnormals, 1);
    if mod(iter,50) == 1
        tangentgradnorm = sqrt(max(sum(tangentgrad .* tangentgrad, 1)))  
    end
% 
    directions = tangentgrad./sqrt(sum(tangentgrad.*tangentgrad,1)); 
    step = sqrt(min(reshape(knn_norms_squared,k_value,[]),[],1));
    cnf_tentative = cnf(:,1:N_moving) +...
                            directions(:,1:N_moving).*step/(offset+iter-1); 
%% Pullback to surface
    h =   surfF(cnf_tentative) ; 
    ngrad = gradF(cnf_tentative);
    while max(abs( h )) > 1e-4    
         cnf_tentative = cnf_tentative - ngrad .* h ./ sum(ngrad .* ngrad, 1);
         h =   surfF(cnf_tentative)  ;
         ngrad = gradF(cnf_tentative);
    end
    cnf = cnf_tentative;
end
toc



% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %                                 
% 
% f_ind = 0;
% %% Triangulate the cnf using surface normals
% for i=1:N
%     nu = null( gradF(cnf(1,i), cnf(2,i), cnf(3,i) )');
%     flatNN = nu \ ( cnf(:,i) - cnf(:,IDX(:,i)) );
%     tri = delaunayTriangulation(flatNN(1,:)',flatNN(2,:)');
%     nbrFaces = tri.vertexAttachments(1);
%     nbrVertsFlat =  tri.ConnectivityList(nbrFaces{:},:);
%     nbrVerts = reshape(IDX(uint8(nbrVertsFlat),i), [],3);
%     faces(f_ind+1:f_ind+size(nbrVerts,1),:) = nbrVerts;
%     f_ind = f_ind+size(nbrVerts,1);
% end
% faces = faces(logical(faces(:,1)),:);
% faces = sort(faces,2);
% faces = unique(faces,'rows');

