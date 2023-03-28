% Mode S1 S2
figure
subplot(1,2,1)
hold on
h3=plot(Ss1(1:end-1), freqs1);
set(h3, 'linewidth', 2)
% h4 =errorbar (xx1, input1, epsilon, ".k");
set(gca, 'fontsize', 14)
ylim([0 200 ] )
yticks([200])
xticks([mean(Ss1)])
xlabel('\it X_1')
pbaspect([1 1 1])
xlim([min(Ss1) max(Ss1) ])
%
subplot(1,2,2)
hold on
h4=plot(Ss2(1:end-1), freqs2);
set(h4, 'linewidth', 2)
set(gca, 'fontsize', 14)
ylim([0 200 ] )
yticks([200])
xticks([mean(Ss2)])
xlabel('\it X_2')
pbaspect([1 1 1])
xlim([min(Ss2) max(Ss2) ])

