function pt_analyzer(cnf, in_domainF)
%PT_ANALYZER
% pt_analyzer(cnf, in_domainF)
% Given a (dim)x(number of points)-array, will determine its separation
% distance, distribution of distances to the nearest neighbors, and radii
% of the largest holes. Provided an indicator function of the restricting
% domain, in_domainF, will detect a subset of nodes, that are at most
% (separation_distance/2) away from the boundary in l1-metric. 
% in_domainF is expected as 
%    in = in_domainF(x, y, z),
% where 'in' is a logical array of the same size as 'x', and 'x', 'y', 'z'
% are the respective coordinates 
%    See also KNNSEARCH, VORONOIN.

path_old = pwd;
path_new = char(mfilename('fullpath'));
cd(path_new(1:end-11))
addpath .
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
if length(size(cnf)) ~= 2
    fprintf( 'Sorry, can only parse a (dim)x(number of points) array.\n')
    return
end
[dim, N] = size(cnf);

format long;
% bins = 200;
% binwidth = .0002;
adjacency = 4*dim+1;
    colors = [[138,186,195]
    [86,54,41]
    [169,191,160]
    [78,55,75]
    [205,181,157]
    [37,63,81]
    [207,175,195]
    [38,66,52]
    [162,169,199]
    [65,63,35]
    [192,146,149]
    [57,110,125]
    [159,136,113]
    [88,96,123]
    [123,144,114]
    [135,118,142]
    [70,101,86]
    [125,88,90]
    [104,146,143]
    [110,99,73]]/256;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
fprintf( 'Analyzing a set of cardinality %d in %d-dimensional space.\n', N, dim);
% % % % % % % % % % SEPARATION OF THE WHOLE NODE SET % % % % % % % % % % %
[~, Dcnf] = knnsearch(cnf', cnf', 'k', adjacency+1);
Dcnf = Dcnf(:,2:end);     % the first column contains only zeros
separation_all = min(Dcnf(:,1));
fprintf( 'The separation distance of this node set:\t%d\n',...
        separation_all);
% % % % % % % % % % SEPARATION OF THE SURFACE NODE SET % % % % % % % % % % 
if exist('in_domainF', 'var') && isa(in_domainF,'function_handle')
    CNF = repmat(cnf,2*dim,1);
    e=[eye(dim) -eye(dim)];
    shifted=bsxfun(@plus, separation_all*e(:)/2,CNF);
    shifted = reshape(shifted,3,[]);
    indices = ~in_domainF( shifted(1,:), shifted(2,:), shifted(3,:));
    indices = reshape(indices,6,[]);
    II=sum(indices,1);
    I=logical(II);
    fprintf( 'The number of nodes l1-close to the surface of the domain:\t%d\n',...
        sum(I));
    cnfsurf=cnf(:,I);          
    [~, Dsurf] = knnsearch(cnfsurf', cnfsurf', 'k', adjacency+1);
    Dsurf = Dsurf(:,2:end);     % the first column contains only zeros
else
     in_domainF = 1;
end
% % % % % % % % % % % % % % % % % % HOLE RADII % % % % % % % % % % % % % %
[V,~] = voronoin(cnf');
if isa(in_domainF,'function_handle')
    V = V(in_domainF(V(:,1),V(:,2),V(:,3)),:);       % Voronoi centers inside the shell
else
    fprintf('No domain indicator function provided, proceeding anyway...\n');
end
[~, holedists] = knnsearch(cnf',V,'k',dim+1);
fprintf('The largest increase from 1st to %d-th deepest hole is:\t %3.6f\n'...
    ,dim+1, max(abs(holedists(:,1)-holedists(:,dim+1))) );
fprintf('(should be zero up to roundoff error, if a domain function is used)\n');


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % PLOTTING % % % % % % % % % % % % % % % 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
if isa(in_domainF,'function_handle')
                                figure(30);
    msize = ceil(max(1, 22-3.5*log10(size(cnfsurf,2)) ));
    plot3(cnfsurf(1,:), cnfsurf(2,:), cnfsurf(3,:),  '.b','MarkerSize',msize);
    axis vis3d;
    pbaspect([1 1 1])
    daspect([1 1 1])
    title('Surface nodes')
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                                figure(50);
pbaspect([1 1 1])
daspect([1 1 1])
hold on;
for i=1:adjacency
    cnf_nearest_neighbors = histogram(Dcnf(:,i));%,bins,'BinWidth',binwidth
    cnf_nearest_neighbors.EdgeAlpha=0;
    cnf_nearest_neighbors.FaceAlpha = .4;
    cnf_nearest_neighbors.Normalization = 'probability';
    cnf_nearest_neighbors.EdgeColor = colors(i,:);
    if isa(in_domainF,'function_handle')
        cnfsurf_nearest_neighbors = histogram(Dsurf(:,i)); %,bins,'BinWidth',binwidth
        cnfsurf_nearest_neighbors.EdgeAlpha=1;   
        cnfsurf_nearest_neighbors.FaceAlpha = .1;
        cnfsurf_nearest_neighbors.Normalization = 'probability';
        cnfsurf_nearest_neighbors.DisplayStyle = 'stairs';  
        cnfsurf_nearest_neighbors.EdgeColor = cnf_nearest_neighbors.EdgeColor;
    end
end
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
deepest_holes = histogram(holedists(:,1)); %,bins,'BinWidth',binwidth
deepest_holes.EdgeAlpha=1;
deepest_holes.DisplayStyle='stairs';
deepest_holes.Normalization = 'probability';
deepest_holes.EdgeColor = 'black';
deepest_holes.LineStyle=':';
deepest_holes.LineWidth=1.5;
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
set(gca,'FontSize',12)
ylabel('Number of nodes','FontSize',24);
xlabel('Distances to the nearest neighbors vs hole radii','FontSize',24);
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% BOXPLOT:
% Input data is specified as a numeric vector or numeric matrix. If x is a
% vector, boxplot plots one box. If x is a matrix, boxplot plots one box 
% for each column of x.
% On each box, the central mark indicates the median, and the bottom and top
% edges of the box indicate the 25th and 75th percentiles, respectively. 
% The whiskers extend to the most extreme data points not considered outliers,
% and the outliers are plotted individually using the '+' symbol.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
                                figure(60);
subplot(1,2,1)
hold on;
boxplot(Dsurf) % ,'PlotStyle','compact'
set(gca,'FontSize',12)
figure(60)
plot(max(Dsurf,[],1));
plot(min(Dsurf,[],1));
leg = legend('Maximal distances for the surface node set','Minimal distances for the surface node set');
leg.FontSize = 16;
leg.Location = 'southeast';
xlim([1 adjacency]);
ylim([0 max(max(Dsurf))]);
%
figure(60)
subplot(1,2,2)
boxplot(Dcnf)
hold on
plot(max(Dcnf,[],1));
plot(min(Dcnf,[],1));
leg = legend('Maximal distances for the whole node set','Minimal distances for the whole node set');
leg.FontSize = 16;
leg.Location = 'southeast';
xlim([1 adjacency]);
ylim([0 max(max(Dsurf))]);
figure(60)
set(gca,'FontSize',12)


% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
cd(path_old)