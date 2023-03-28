% 2022-01-07
figure
subplot(1,2,1)
hold on
for ii = 1:length(ROI_array)
  clear y_now;
  ROI_now = ROI_array(ii,1):  ROI_array(ii,2);
  data_now = input1(ROI_now);
  in_now = ROI_now';
  [tau1, w1, yint] = DataLinearModelZ (data_now, epsilon0);

 y_now = data_now  - tau1(2)*(in_now-in_now(1)+1) - tau1(1)+input1(1);
 h3 = plot(in_now, y_now, '-y')
 h1=errorbar (in_now, y_now, epsilon(ROI_now), ".b");
yynew1(ROI_now) = y_now;
  end
set(gca, 'fontsize', 14)
ylim([min(yynew1)-epsilon0 max(yynew1)+epsilon0 ] )
xlim([1 200])
%
subplot(1,2,2)
hold on
for ii = 1:length(ROI_array)
  clear yp;
  ROI_now = ROI_array(ii,1):  ROI_array(ii,2);
  data_now = input2(ROI_now);
  in_now = ROI_now';
  [tau1, w1, yint] = DataLinearModelZ (data_now, epsilon0);

  y_now = data_now  - tau1(2)*(in_now-in_now(1)+1) - tau1(1)+input2(1);
 h3 = plot(in_now, y_now, '-y')
 h1=errorbar (in_now, y_now, epsilon(ROI_now), ".b");
yynew2(ROI_now) = y_now;
  end
set(gca, 'fontsize', 14)
ylim([min(yynew2)-epsilon0 max(yynew2)+epsilon0 ] )
xlim([1 200])

