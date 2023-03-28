%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-10-10
% 2022-10-11
figure
subplot(4,1,1)
hold on
h1 = plot(R21est, JK_arrayP, '-k', "linewidth", 1)
##h1 = plot(R21est(ind_out_L), JK_array(ind_out_L), '-k', "linewidth", 1)
##h2 = plot(R21est(ind_out_R), JK_array(ind_out_R), '-k', "linewidth", 1)
h2 = plot(R21est(max_T123_ind_more), JK_arrayP(max_T123_ind_more), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
hh = axis
##xx = [R21est(min(max_JK_array_ind)) R21est(min(max_JK_array_ind))]
##yy = [hh(3) JK_arrayP(min(max_JK_array_ind))]
##plot(xx, yy, '--r')
##xx = [R21est(max(max_JK_array_ind)) R21est(max(max_JK_array_ind))]
##yy = [hh(3) JK_arrayP(max(max_JK_array_ind))]
##plot(xx, yy, '--r')
##xticks([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind)) ])
##xticklabels([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind))])
ylim([min(JK_arrayP) max(JK_arrayP)])
xlim([min(R21est) max(R21est)])
ylabel('\it JK')
%
##xx = [R21est(max_JK_array_ind) R21est(max_JK_array_ind)]
##yy = [hh(3) JK_arrayP(max_JK_array_ind)]
##plot(xx, yy, '--r')
xticks([R21est(max_JK_array_ind)  ])
xticklabels([R21est(max_JK_array_ind) ])
xticklabels([])
yticks([min(JK_arrayP) max(JK_arrayP)])
yticklabels([min(JK_arrayP) max(JK_arrayP)])
%
subplot(4,1,2)
hold on
h1 = plot(R21est, Kinv, '-k', "linewidth", 1)
h2 = plot(R21est(max_T123_ind_more), Kinv(max_T123_ind_more), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
ylim([min(Kinv) max(Kinv)])
xlim([min(R21est) max(R21est)])
xticklabels([])
yticklabels([])
yticks([min(Kinv) max(Kinv)])
yticklabels([min(Kinv) max(Kinv)])
ylim([[min(Kinv) max(Kinv)]])
ylabel('\it Kinv_{O}')
hh = axis
##xx = [R21est(min(max_JK_array_ind)) R21est(min(max_JK_array_ind))]
##yy = [hh(3) Kinv(min(max_JK_array_ind))]
##plot(xx, yy, '--r')
##xx = [R21est(max(max_JK_array_ind)) R21est(max(max_JK_array_ind))]
##yy = [hh(3) Kinv(max(max_JK_array_ind))]
##plot(xx, yy, '--r')
xticks([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind)) ])
xticklabels([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind))])
xticklabels([])
%
subplot(4,1,3)
hold on
h1 = plot(R21est, max_mu_array01, '-k', "linewidth", 1)
h2 = plot(R21est(max_T123_ind_more), max_mu_array01(max_T123_ind_more), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
ylim([min(max_mu_array01) max(max_mu_array01)])
xlim([min(R21est) max(R21est)])
xticklabels([])
yticklabels([])
xticks([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind)) ])
xticklabels([R21est(min(max_JK_array_ind)) R21est(max(max_JK_array_ind))])
xticklabels([])
yticks([min(max_mu_array01) max(max_mu_array01)])
yticklabels([min(max_mu_array01) max(max_mu_array01)])
%xlabel('\it k_{12}')
ylabel('\it max \mu')
set(gca, 'linewidth', 1);
set(gca, 'FontSize', 14);
hh = axis
%

alpha_set_mu = find (m_array > alpha_lev_mu*max(m_array) )
R21est(min(alpha_set_mu))
R21est(max(alpha_set_mu))

xxx =  [ R21est(min(alpha_set_mu)) R21est(min(alpha_set_mu)) ]
yyy = [ 0 alpha_lev_mu*max(m_array)  ]
h1 = plot(xxx, yyy, '--r')
%
xxx =  [ R21est(max(alpha_set_mu)) R21est(max(alpha_set_mu))]
yyy = [ 0 alpha_lev_mu*max(m_array) ]
h1 = plot(xxx, yyy, '--r')
##xx = [R21est(min(max_JK_array_ind)) R21est(min(max_JK_array_ind))]
##yy = [hh(3) max_mu_array01(min(max_JK_array_ind))]
##plot(xx, yy, '--r')
##xx = [R21est(max(max_JK_array_ind)) R21est(max(max_JK_array_ind))]
##yy = [hh(3) max_mu_array01(max(max_JK_array_ind))]
##plot(xx, yy, '--r')

subplot(4,1,4)
hold on
h1 = plot(R21est, T123, '-k', "linewidth", 1)
h2 = plot(R21est(max_T123_ind_more), T123(max_T123_ind_more), '-r', "linewidth", 1)
box('off')
ylim([min(T123) max(T123)])
xlim([min(R21est) max(R21est)])
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
xticklabels([])
yticklabels([])
yticks([min(T123) max(T123)])
yticklabels([min(T123) max(T123)])
xlabel('\it R')
ylabel('\it T')
set(gca, 'linewidth', 1);
set(gca, 'FontSize', 14);
hh = axis
##xx = [R21est(min(max_JK_array_ind)) R21est(min(max_JK_array_ind))]
##yy = [hh(3) T123(min(max_JK_array_ind))]
##plot(xx, yy, '--r')
##xx = [R21est(max(max_JK_array_ind)) R21est(max(max_JK_array_ind))]
##yy = [hh(3) T123(max(max_JK_array_ind))]
##plot(xx, yy, '--r')
##xx = [R21est(max_JK_array_ind) R21est(max_JK_array_ind)]
##yy = [hh(3) T123(max_JK_array_ind)]
##plot(xx, yy, '--r')
##xticks([R21est(min(max_JK_array_ind)) R21est(max_JK_array_ind) R21est(max(max_JK_array_ind)) ])
##xticklabels([R21est(min(max_JK_array_ind)) R21est(max_JK_array_ind) R21est(max(max_JK_array_ind))])

xx = [R21est(min(max_T123_ind_more)) R21est(min(max_T123_ind_more))]
yy = [hh(3) T123(min(max_T123_ind_more))]
plot(xx, yy, '--r')
xx = [R21est(max(max_T123_ind_more)) R21est(max(max_T123_ind_more))]
yy = [hh(3) T123(max(max_T123_ind_more))]
plot(xx, yy, '--r')
xticklabels([])
xticks([])

% 2023-01-19

%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
