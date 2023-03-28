%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-10-10
% 2022-10-11
% 2022-12-26
figure
subplot(3,1,1)
hold on
h1 = plot(R21est, JK_array, '-k', "linewidth", 1)
##h1 = plot(R21est(ind_out_L), JK_array(ind_out_L), '-k', "linewidth", 1)
##h2 = plot(R21est(ind_out_R), JK_array(ind_out_R), '-k', "linewidth", 1)
h2 = plot(R21est(max_JK_array_ind ), JK_array(max_JK_array_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
hh = axis
xx = [R21est(min(max_JK_array_ind)) R21est(min(max_JK_array_ind))]
yy = [hh(3) JK_array(min(max_JK_array_ind))]
plot(xx, yy, '--r')
xticks([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind)) ])
xticklabels([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind))])
ylim([min(JK_array) max(JK_array)])
xlim([min(R21est) max(R21est)])
yticks([min(JK_array) max(JK_array)])
yticklabels([min(JK_array) max(JK_array)])
ylabel('\it JK')
%
xticks([R21est(max_JK_array_ind)  ])
xticklabels([R21est(max_JK_array_ind) ])
ylim([min(JK_array) max(JK_array)])
%
subplot(3,1,2)
hold on
h1 = plot(R21est, k_array, '-k', "linewidth", 1)
h2 = plot(R21est(max_JK_array_ind ), k_array(max_JK_array_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
xticklabels([])
yticklabels([])
xlim([min(R21est) max(R21est)])
yticks([min(k_array) max(k_array)])
yticklabels([min(k_array) max(k_array)])
ylim([min(k_array) max(k_array)])
ylabel('\it K_{O}')
hh = axis
xx = [R21est(min(max_JK_array_ind)) R21est(min(max_JK_array_ind))]
yy = [hh(3) k_array(min(max_JK_array_ind))]
plot(xx, yy, '--r')
xticks([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind)) ])
xticklabels([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind))])
%
subplot(3,1,3)
hold on
h1 = plot(R21est, m_array, '-k', "linewidth", 1)
h2 = plot(R21est(max_JK_array_ind ), m_array(max_JK_array_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
xlim([min(R21est) max(R21est)])
xticklabels([])
yticklabels([])
xticks([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind)) ])
xticklabels([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind))])
yticks([min(m_array) max(m_array) ])
yticklabels([min(m_array) max(m_array) ])
ylim([min(m_array) max(m_array) ])
xlabel('\it R')
ylabel('\it max \mu')
set(gca, 'linewidth', 1);
set(gca, 'FontSize', 14);
hh = axis
xx = [R21est(min(max_JK_array_ind)) R21est(min(max_JK_array_ind))]
yy = [hh(3) m_array(min(max_JK_array_ind))]
plot(xx, yy, '--r')

% 2023-01-19
alpha_lev = 0.985 %0.999%
alpha_set = find (m_array > alpha_lev*max(m_array) )
R21est(min(alpha_set))
R21est(max(alpha_set))
%
xxx =  [ R21est(min(alpha_set)) R21est(min(alpha_set)) ]
yyy = [ 0 alpha_lev*max(m_array)  ]
h1 = plot(xxx, yyy, '--r')
%
xxx =  [ R21est(max(alpha_set)) R21est(max(alpha_set)) ]
yyy = [ 0 alpha_lev*max(m_array) ]
h1 = plot(xxx, yyy, '--r')
%

%figure_name_out=strcat('PgammaPhTotalJKOskorbinMu', '.png')
%figure_name_out=strcat(FNstr,'JKOskorbinMu', '.png')
%print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
