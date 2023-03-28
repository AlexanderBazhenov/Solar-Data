% 2022-01-07
figure
subplot(1,2,1)
hold on
for ii = 1:length(ROI_array)
  clear yp;
  ROI_now = ROI_array(ii,1):  ROI_array(ii,2);
  data_now = input1(ROI_now);
  in_now = ROI_now;
  [tau1, w1, yint] = DataLinearModelZ (data_now, epsilon0);
  h1=errorbar (in_now, data_now, epsilon(ROI_now), ".b");
%
  yp= tau1(1) +tau1(2)*(in_now-in_now(1)+1);


h2=plot(in_now, yp, '-y')

  end
set(gca, 'fontsize', 14)
ylim([min(input1)-epsilon0 max(input1)+epsilon0 ] )
xlim([1 200])
subplot(1,2,2)
hold on
for ii = 1:length(ROI_array)
  clear yp;
  ROI_now = ROI_array(ii,1):  ROI_array(ii,2);
  data_now = input2(ROI_now);
  in_now = ROI_now;
  [tau1, w1, yint] = DataLinearModelZ (data_now, epsilon0);
  h1=errorbar (in_now, data_now, epsilon(ROI_now), ".b");
%
  yp= tau1(1) +tau1(2)*(in_now-in_now(1)+1);


h2=plot(in_now, yp, '-y')

  end
set(gca, 'fontsize', 14)
ylim([min(input2)-epsilon0 max(input2)+epsilon0 ] )
xlim([1 200])

