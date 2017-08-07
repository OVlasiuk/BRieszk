function surf_voronoi(cnf, surfF, gradF)
%SURF_VORONOI
% Approximate Voronoi diagram on the level surface.
% cnf -- 3x(num_pts), the node set to be processed
% surfF -- the surface is defined by surfF(x,y,z) == 0;
% gradF -- a function handle to evaluate the gradient to the surface at
%   a point (x,y,z);

k_value = 100;
N = size(cnf,2);

[IDX, ~] = knnsearch(cnf', cnf', 'k', k_value+1);
IDX = IDX';

faces = zeros(10 * k_value * N, 3);         % not a real estimate, of course; 
                                            % a bad upper bound is 
                                            % n choose 3
C_Tri = round(rand(k_value * N, 1)*3);
f_ind = 0;
for i=1:N
    nu = null( gradF(cnf(1,i), cnf(2,i), cnf(3,i) )');
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

C_Tri = C_Tri(logical(faces(:,1)),:);
T = triangulation(faces, cnf(1,:)', cnf(2,:)', cnf(3,:)');
% % Plot the triangulation
% figure
% TS = trisurf(T);
% % TS = trimesh(T);
% % set(gca,'CLim',[min(C_Tri), max(C_Tri)]);
% % set(TS,'FaceColor','flat',...
% %        'FaceVertexCData',C_Tri,...
% %        'CDataMapping','scaled',...
% %        'FaceAlpha', .2,...
% %        'EdgeAlpha',1);
% set(TS,'FaceColor','flat',...
%        'FaceAlpha', 1,...
%        'EdgeAlpha',1);
% pbaspect([1 1 1])
% daspect([1 1 1])
% set(gca, 'Clipping', 'off')
% axis vis3d

% FaceVertexAlphaData!
 
voronois = T.circumcenter;         % long thin
vorDiagramVerts = [cnf'; voronois];
vorDiagramFaces = zeros(3*size(faces,1), 3);
C_TriVor = zeros(3*size(faces,1),1);
fill_i = 0;
for i=1:N
    adjFaces = T.vertexAttachments{i};
    nu = null( gradF(cnf(1,i), cnf(2,i), cnf(3,i) )');
%   Order the adjacent faces by angles between the centroid and the origin.
    adjVerts = T.ConnectivityList(adjFaces,:);
    adjVertsCoords = cnf(:, adjVerts') - cnf(:, i);
    adjVertsCoordsProjected = nu \ adjVertsCoords;
    faceGroupedVerts = reshape(adjVertsCoordsProjected, 2, [], size(adjVerts,1));
    faceCentroids = squeeze(sum(faceGroupedVerts, 2));
    [~, sortInd] = sort( atan2(faceCentroids(1,:), faceCentroids(2,:)) );
    adjFaces = adjFaces(sortInd);
%     
    vorDiagramFaces( fill_i+1:fill_i+numel(adjFaces), :) =...
        [N+adjFaces; circshift(N+adjFaces,-1); i*ones(size(adjFaces))]';
    C_TriVor( fill_i+1:fill_i+numel(adjFaces) ) =...
        numel(adjFaces) * ones(numel(adjFaces),1);
    fill_i = fill_i + numel(adjFaces);
end
% vorEdges = zeros(fill_i/2, 2);


    
TVor = triangulation(vorDiagramFaces, vorDiagramVerts(:,1), vorDiagramVerts(:,2), vorDiagramVerts(:,3));
f = figure
TSVor = trisurf(TVor);
set(gca,'CLim',[min(C_TriVor), max(C_TriVor)]);
set(TSVor,'FaceColor','flat',...
       'FaceVertexCData',C_TriVor,...
       'CDataMapping','scaled',...
       'EdgeAlpha',0);
pbaspect([1 1 1])
daspect([1 1 1])
set(gca, 'Clipping', 'off')
axis vis3d
camzoom(1.2)
f.PaperType = 'a3';
print(f, 'voronoi','-djpeg','-r600')

% hold on;

% TMesh = trimesh(TVor);
% set(TMesh, 'EdgeColor','k','FaceAlpha',0);

% 
% % f = @(x) x(1).^2 .*(x(1).^2 - 5) + x(2).^2 .*(x(2).^2 - 5) + x(3).^2 .*(x(3).^2 - 5) + 11;
% % N = 10000;
% % x=zeros(3,N);
% % x0 = 5*rand(3,N)-2.5;
% % for i=1:N
% % x(:,i) = fsolve(f, x0(:,i),optimoptions('fsolve','Display','off'));
% % end
% % plot3(x(1,:),x(2,:),x(3,:),'.b')
% % pbaspect([1 1 1])
% % daspect([1 1 1])
% % axis vis3d
% 
% 
% % g = @(x,y,z) (1 - z.^2) .* (2* y .* (y.^2 - 3 *x.^2) - 9* z.^2 + 1) + (x.^2 + y.^2).^2;
% % fs = fimplicit3(g)
% % axis vis3d
% % g = @(x,y,z) x.^2+y.^2+z.^2+15 -8*sqrt(x.^2+y.^2) % torus 

% % f = @(x,y,z) x.^2 .*(x.^2 - 5) + y.^2 .*(y.^2 - 5) + z.^2 .*(z.^2 - 5) + 11;
% % fs = fimplicit3(f)
% % axis vis3d

