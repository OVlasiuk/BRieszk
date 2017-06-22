function ratio = dcompare(pts,densityF, plotit)
%DCOMPARE
% function dcompare(pts,densityF)
% Display statistics about how the radial density (aka distance to the nearest neighbor)
% compares with the values of function densityF at the points 'pts'. Must be vectors in 3-d space.
%
% pts -- 3x(numpts) array
% densityF -- function handle
if size(pts,1) ~= 3
    pts = pts';
end

[~, D] = knnsearch(pts', pts', 'k', 2);
rdens_cnf = D(:,2);
rdens_fun = densityF(pts);
ratio = rdens_fun./rdens_cnf';
diff = abs(rdens_fun - rdens_cnf');
radii = sqrt( sum( pts .*pts,1 ) );
maxdiff = max(diff);
meandiff = mean(diff);
minratio = min(ratio);
maxratio = max(ratio);
quantile5 = quantile(ratio,0.05);
quantile95 = quantile(ratio,0.95);
meanratio = mean(ratio);
varratio = var(ratio);
fprintf('\nmaxdiff\t\tmeandiff\n');
fprintf('%3.6f\t%3.6f\n\n',maxdiff,meandiff)
fprintf('minratio\tmaxratio\tquantile5\tquantile95\n');
fprintf('%3.6f\t%3.6f\t%3.6f\t%3.6f\t\n\n', minratio,maxratio,quantile5,quantile95)
fprintf('meanratio\tvarratio\n')
fprintf('%3.6f\t%3.6f\n',meanratio,varratio)

if exist('plotit','var') && (plotit=='y' || plotit=='Y' || plotit==1)
    msize = ceil(max(1, 22-5*log10(size(pts,2)) ));
    figure(3);
    plot(radii,ratio,'.k', 'MarkerSize',4)
    hold on;
    plot(radii,diff,'.g', 'MarkerSize',4)
    set(gca,'FontSize',12)
    xlabel('Radius {\bf\it{N}}','FontSize',24);
    ylabel('\rho({\bf\it{N}})/\Delta({\bf\it{N}})','FontSize',24);
    
    
    figure(4)   
    plot3(pts(1,ratio>quantile95),pts(2,ratio>quantile95),pts(3,ratio>quantile95),'.k','MarkerSize',msize)
    hold on;
    plot3(pts(1,ratio<quantile5),pts(2,ratio<quantile5),pts(3,ratio<quantile5),'.r','MarkerSize',msize)
    axis vis3d;
    daspect([1 1 1]);
    pbaspect([1 1 1]);
end
