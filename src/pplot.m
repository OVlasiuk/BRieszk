function pplot(cnf)
if size(cnf,1) > size(cnf,2)
    cnf = cnf';
end
msize = ceil(max(1, 22-3.5*log10(size(cnf,2)) ));
colormap white
plot3(cnf(1,:),cnf(2,:),cnf(3,:),'.k','MarkerSize',msize)
pbaspect([1 1 1])
daspect([1 1 1])
set(gca, 'Clipping', 'off')
axis vis3d
box on
grid on


%     scatter3(cnf(1,:),cnf(2,:),cnf(3,:),'MarkerFaceColor','r',...
%         'MarkerFaceAlpha',.2,'MarkerEdgeAlpha',.2); 
