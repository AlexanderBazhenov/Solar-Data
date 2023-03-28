close all
figure
subplot(1,2,1)
h3=errorbar (xx1, input1, epsilon, ".b");
set(gca, 'fontsize', 14)
ylim([min(input1)-epsilon0 max(input1)+epsilon0 ] )
xlim([0 length(xx1)+1])
xlabel('\it n')
pbaspect([1 2 1])
%
subplot(1,2,2)
h4=errorbar (xx1, input2, epsilon, ".r");
set(gca, 'fontsize', 14)
ylim([min(input2)-epsilon0 max(input2)+epsilon0 ] )
xlim([0 length(xx1)+1])
xlabel('\it n')
pbaspect([1 2 1])
##figure_name_out=strcat(FNstr,'Inte','2Axis', '.png')
##print('-dpng', '-r300', figure_name_out), pwd




