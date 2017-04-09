function riesz_sphere

k_value = 20;  
dim = 3;
s = 6;
N = 10000;
repel_steps = 200;
% cnf = randn(dim,N); 
% cnf = cnf./sqrt(sum(cnf.^2,1));
% - this is faster than normc in Matlab 2016b
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
close all;
pbaspect([1 1 1])
colormap(winter)
[x,y,z] = sphere(30);
mesh(x,y,z)
hold on
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',8)


for cycle=1:10
    [IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
    IDX = IDX(:,2:end)';             % drop the trivial first column in IDX
    step = min(D(:,2));
    tic
    for iter=1:repel_steps
        cnf_repeated = reshape(repmat(cnf,k_value,1),dim,[]);
        knn_differences = cnf_repeated - cnf(:,IDX);
        knn_norms = sqrt(sum(knn_differences.^2,1));
        riesz_weights = knn_norms.^(-s-1);
        directions = bsxfun(@times,riesz_weights,knn_differences);
        directions = reshape(directions, dim, k_value, []);
        normals = reshape(sum(directions,2), dim, []);
%         normals = directions./sqrt(sum(directions.^2,1));
        tangents = normals-bsxfun(@times,sum(cnf.*normals,1),cnf); 
        % computing scalar products only works because we are on the unit
        % sphere
        tangents = tangents/max(sqrt(sum(tangents.^2,1)));
        cnf = cnf + tangents * step/3/cycle;%/iter
        cnf = cnf./sqrt(sum(cnf.^2,1));
    end
    toc
end
[IDX, D] = knnsearch(cnf', cnf', 'k', k_value+1);
step = min(D(:,2));
plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',7)
% dlmwrite('cnf.out',cnf','delimiter','\t');