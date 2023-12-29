function Image = Plot(pr_, edge_, Image, dp, pos_mask, neg_mask)
% Plot
% Inputs:
%   - Image: Image magnetogram.
%   - pr_: Page-Rank.
%   - edge_: Degree of Node.
%   - dp: The matrix includes the total degree of each patch, x, and y of the center of the patches.
%   - pos_mask and neg_mask: boundaries.

subplot(1,3,3),imshow (pr_,[-60 60])
imagesc(pr_)
axis xy
colorbar
colormap(pink(256));
xlabel('Solar-X','FontSize',22,'FontName','Times New Roman')
set(gca,'FontName','Times New Roman','FontSize',22)
set(gca,'FontName','Times New Roman','FontSize',22);
title('Page-Rank Map')

subplot(1,3,1),imshow(Image,[-60 60])
axis on
axis xy
ylabel('Solar-Y','FontSize',22,'FontName','Times New Roman')
xlabel('Solar-X','FontSize',22,'FontName','Times New Roman')
title('HMI Image')
set(gca,'FontName','Times New Romhisan','FontSize',22)
set(gca,'FontName','Times New Roman','FontSize',22);


subplot(1,3,2),imshow(edge_,[-60 60])
title('Degree of Node Map')
set(gca,'FontName','Times New Roman','FontSize',22)
set(gca,'FontName','Times New Roman','FontSize',22);
axis on
axis xy
hold on
contour(pos_mask,1,'LineWidth',2,'LineColor','r')
contour(neg_mask,1,'LineWidth',2,'LineColor','r')
set(gca, 'FontName', 'Times New Romhisan', 'FontSize', 22)
set(gca, 'FontName', 'Times New Roman', 'FontSize', 22);
xlabel('Solar-X', 'FontSize', 22, 'FontName', 'Times New Roman')
title('Degree of node Map')
for i=1:length(dp(:,1))
    text(round(dp(i,2)), round(dp(i,3)), sprintf('%d', i), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle','color','b',...
        'FontSize',14,'FontName','Times New Roman',...
        'FontWeight','bold');
    hold on
end


end
