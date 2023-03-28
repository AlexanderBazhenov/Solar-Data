figure
subplot(1,2,1)
hold on
h3=errorbar (xx1',  yynew1, epsilon, ".b");
% h4 =errorbar (xx1, input1, epsilon, ".k");
set(gca, 'fontsize', 14)
ylim([min(yynew1)-epsilon0 max(yynew1)+epsilon0 ] )
xlim([0 length(xx1)+1])
xp = [0 length(xx1)+1]
yp =[ mid(inner1) mid(inner1)]
plot( xp, yp, '-k')
box('on')
xlabel('\it n')
pbaspect([1 2 1])
%
subplot(1,2,2)
hold on
h4=errorbar (xx1', yynew2, epsilon, ".r");
set(gca, 'fontsize', 14)
ylim([min(yynew2)-epsilon0 max(yynew2)+epsilon0 ] )
xlim([0 length(xx1)+1])
xlabel('\it n')
xp = [0 length(xx1)+1]
yp =[ mid(inner2) mid(inner2)]
box('on')
plot( xp, yp, '-k')
pbaspect([1 2 1])
figure_name_out=strcat(FNstr,'InteConst','2Axis', '.png')
print('-dpng', '-r300', figure_name_out), pwd
