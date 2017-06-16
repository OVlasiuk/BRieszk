function dcompare(pts,densityF)
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
