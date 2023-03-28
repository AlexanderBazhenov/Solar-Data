% 2022-12-30
%close all
[tau1, w1, yint] = DataLinearModel (input1, epsilon0);
[tauz, wz, yinz] = DataLinearModelZ (input1, epsilon0);
figure
subplot(1,2,1)
h3=errorbar (xx1, input1, epsilon, ".b");
set(gca, 'fontsize', 14)
ylim([min(input1)-epsilon0 max(input1)+epsilon0 ] )
%
yp=tau1(1)+tau1(2)*xx1;
hold on
h2=plot(xx1, yp, '-r')
%
yp=tauz(1)+tauz(2)*xx1;
hold on
h2=plot(xx1, yp, '-y')
%
xlim([0 length(xx1)+1])
xlabel('\it n')
pbaspect([1 2 1])
%
[tau2, w2, yint] = DataLinearModel (input2, epsilon0);
[tauz2, wz2, yinz2] = DataLinearModelZ (input2, epsilon0);
subplot(1,2,2)
h4=errorbar (xx1, input2, epsilon, ".r");
set(gca, 'fontsize', 14)
ylim([min(input2)-epsilon0 max(input2)+epsilon0 ] )
xlim([0 length(xx1)+1])
%
yp=tau2(1)+tau2(2)*xx1;
hold on
h2=plot(xx1, yp, '-b')
%
yp=tauz2(1)+tauz2(2)*xx1;
hold on
h2=plot(xx1, yp, '-g')
%
xlabel('\it n')
pbaspect([1 2 1])
##figure_name_out=strcat(FNstr,'Inte','2Axis','LinMode', '.png')
##print('-dpng', '-r300', figure_name_out), pwd




