function [vorFig, triFig] = f_vorsurf(cnf, gradF, densityF, varargin)
%SURF_VORONOI
% [vorFig, triFig] = surf_voronoi(cnf, gradF)
% Approximate Voronoi diagram on the level surface.
% INPUT:
% cnf -- 3x(num_pts), the node set to be processed
% gradF -- a function handle to evaluate the gradient to the surface at
%   a point (x,y,z);
% gradF -- a function handle for the density of the distribution, will be
%           used for coloring the triangulation;
% Optional argument name/value pairs:
% Name          Value
%
% 'k'           number of nearest neighbors used; default: 25
% 'defaults'    pass 'true' to use the default gradient/density;
%               default: false

% OUTPUT:
% vorFig -- handle to the figure with the surface Voronoi diagram;
% triFig -- handle to the surface triangulation with vertices at cnf.
% Both are returned with the .Visible attribute set to 'off'.
%
% 
% EXAMPLE of a configuration on a complicated surface:
% warning('off','optim:fsolve:NonSquareSystem')
% f = @(x) x(1).^2 .*(x(1).^2 - 5) + x(2).^2 .*(x(2).^2 - 5) +...
% x(3).^2 .*(x(3).^2 - 5) + 11;
% N = 1000;
% x=f_cnfinit(N, f);
% plot3(x(1,:),x(2,:),x(3,:),'.b')
% pbaspect([1 1 1])
% daspect([1 1 1])
% axis vis3d
% 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

pnames = { 'k', 'defaults'};
dflts =  { 25 false };
[k_value, defaults, ~] =...
     internal.stats.parseArgs(pnames, dflts, varargin{:});
if defaults
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
end
N = size(cnf,2);
msize = ceil(max(1, 22-6*log10(size(cnf,2)) ));

[IDX, ~] = knnsearch(cnf', cnf', 'k', k_value+1);
IDX = IDX';

faces = zeros(10 * k_value * N, 3);         % not a real estimate, of course; 
                                            % a bad upper bound is 
                                            % n choose 3                         
f_ind = 0;
%% Triangulate the cnf using surface normals
for i=1:N
    nu = null( gradF(cnf(:,i))');
    flatNN = nu \ ( cnf(:,i) - cnf(:,IDX(:,i)) );
    tri = delaunayTriangulation(flatNN(1,:)',flatNN(2,:)');
    nbrFaces = tri.vertexAttachments(1);
    nbrVertsFlat =  tri.ConnectivityList(nbrFaces{:},:);
    nbrVerts = reshape(IDX(uint8(nbrVertsFlat),i), [],3);
    faces(f_ind+1:f_ind+size(nbrVerts,1),:) = nbrVerts;
    f_ind = f_ind+size(nbrVerts,1);
end
faces = faces(logical(faces(:,1)),:);
faces = sort(faces,2);
faces = unique(faces,'rows');

%% Draw the triangulation
T = triangulation(faces, cnf(1,:)', cnf(2,:)', cnf(3,:)');
if isa(densityF,'function_handle')
    C_Tri = densityF( cnf(:,faces(:)) );
    C_Tri = reshape(C_Tri,[],3);
    C_Tri = mean(C_Tri,2);
else
    C_Tri = round(rand(size(faces,1), 1) * 4);
end

%% Plot the triangulation
triFig = figure;
set(triFig, 'Visible', 'off');
TS = trisurf(T);
set(gca,'CLim',[min(C_Tri), max(C_Tri)]);
set(TS,'FaceColor','flat',...
       'FaceVertexCData',C_Tri,...
       'CDataMapping','scaled',...
       'FaceAlpha', 1,...
       'EdgeAlpha',1);
% set(TS,'FaceColor','flat',...
%        'FaceAlpha', 1,...
%        'EdgeAlpha',1);
pbaspect([1 1 1])
daspect([1 1 1])
% set(gca, 'Clipping', 'off')
axis vis3d

%% Determine faces of the Voronoi diagram 
voronois = T.circumcenter;         % long thin
vorDiagramVerts = [cnf'; voronois];
vorDiagramFaces = zeros(3*size(faces,1), 3);
vorDiagramEdges = zeros(3*size(faces,1), 2);
C_TriVor = zeros(3*size(faces,1),1);
fill_i = 0;
for i=1:N
    adjFaces = T.vertexAttachments{i};
    nu = null( gradF(cnf(:,i) )');
%   Order the adjacent faces by angles between the centroid and the origin.
    adjVerts = T.ConnectivityList(adjFaces,:);
    adjVertsCoords = cnf(:, adjVerts') - cnf(:, i);
    adjVertsCoordsProjected = nu \ adjVertsCoords;
    faceGroupedVerts = reshape(adjVertsCoordsProjected, 2, [], size(adjVerts,1));
    faceCentroids = sum(faceGroupedVerts, 2);
    [~, sortInd] = sort( atan2(faceCentroids(1,1,:), faceCentroids(2,1,:)) );
    adjFaces = adjFaces(sortInd);
%     
    vorDiagramFaces( fill_i+1:fill_i+numel(adjFaces), :) =...
        [N+adjFaces; circshift(N+adjFaces,-1); i*ones(size(adjFaces))]';
    vorDiagramEdges( fill_i+1:fill_i+numel(adjFaces), :) =...
        [N+adjFaces; circshift(N+adjFaces,-1); ]';
    C_TriVor( fill_i+1:fill_i+numel(adjFaces) ) =...
        numel(adjFaces) * ones(numel(adjFaces),1);
    fill_i = fill_i + numel(adjFaces);
end


%% Draw the Voronoi diagram
TVor = triangulation(vorDiagramFaces, vorDiagramVerts(:,1), vorDiagramVerts(:,2), vorDiagramVerts(:,3));
vorFig = figure;
set(vorFig, 'Visible', 'off');
TSVor = trisurf(TVor);
set(gca,'CLim',[4, 7.5]);
set(TSVor,'FaceColor','flat',...
       'FaceVertexCData',C_TriVor,...
       'CDataMapping','scaled',...
       'EdgeAlpha',0);
hold on;
TEdges = trisurf(vorDiagramEdges(:,[1 1 2]),vorDiagramVerts(:,1), vorDiagramVerts(:,2), vorDiagramVerts(:,3));
set(TEdges,'EdgeAlpha',.3 ...      % determined the edge transparency 
       );
plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',msize)

pbaspect([1 1 1])
daspect([1 1 1])
view(3)
% campos([-45.6343  -59.4699   43.2793])
% set(gca, 'Clipping', 'off')
% set(gca,'xtick','')
% set(gca,'ytick','')
% set(gca,'ztick','')
% axis off
axis vis3d
% camzoom(1.9);


% % g = @(x,y,z) (1 - z.^2) .* (2* y .* (y.^2 - 3 *x.^2) - 9* z.^2 + 1) + (x.^2 + y.^2).^2;
% % fs = fimplicit3(g)
% % axis vis3d
% % g = @(x,y,z) x.^2+y.^2+z.^2+15 -8*sqrt(x.^2+y.^2) % torus 

% % f = @(x,y,z) x.^2 .*(x.^2 - 5) + y.^2 .*(y.^2 - 5) + z.^2 .*(z.^2 - 5) + 11;
% % fs = fimplicit3(f)
% % axis vis3d

% for i=10:10:3650
% cnf = dlmread([int2str(i) 'it.txt'])';
% f = surf_voronoi(cnf,g,grad);
% print(f, [int2str(i) 'it'],'-djpeg','-r800')
% close
% end
