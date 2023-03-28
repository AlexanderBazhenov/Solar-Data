% 2022-10-11
% PgammaData.m
%
clear all
close all

dirroot = 'D:\Data\ST\2022\T\'
%
dirroot ='e:\Users\Public\Documents\ST\2022\T\'
% 2022-06-16
% ki
cd(dirroot), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Phys Lett 1992
Pgamma1992Data;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

OxfordBlue = [0, .33, .71]
RoyalMail = 4.58*[0.218, .032, 0.042]
Pantone = 3*[0.128, 0.140, 0.036]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pkg load interval
PgammaStd = midrad(PgammaPh, PgammaPhStat)
PgammaComptonStd =midrad(PgammaCompton, PgammaComptonStat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JK_array = []
max_mu_array = []
k_array = []
%
kcorrCompt_array=1:0.025:2
%
for jj = 1:length(kcorrCompt_array)
  kc_now=kcorrCompt_array(jj)
  X = [ PgammaStd; kc_now*PgammaComptonStd];
  [JK_now, max_mu_now, k_now] = SetMeasures(X);
  JK_array = [JK_array JK_now];
  max_mu_array = [max_mu_array max_mu_now];
  k_array = [k_array k_now];
end
% save optimization arrays
save PgammaArrays.mat JK_array max_mu_array k_array kcorrCompt_array
% load PgammaArrays
% good optimization interval
max_mu_ind = find (max_mu_array >= max(max_mu_array))
Oskorbin_ind = find (k_array <= min(k_array)+0.0001)
max_mu_Oskorbin_ind = intersect(max_mu_ind, Oskorbin_ind)
ind_out = setdiff(1:length(kcorrCompt_array), max_mu_Oskorbin_ind )
% out of good
[indmax iinmax2] = max(diff(ind_out))
ind_out_L = ind_out(1:iinmax2)
ind_out_R = setdiff( ind_out, ind_out_L )
%
kcorrCompt_array(max_mu_Oskorbin_ind )
%
[max_JK_array, max_JK_array_ind] = max(JK_array)
[min_JK_array, min_JK_array_ind] = min(JK_array)
k_array_opt =kcorrCompt_array(max_JK_array_ind)
max_mu_array_opt=max_mu_array(max_JK_array_ind)


%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-10-10
% 2022-10-11
figure
subplot(3,1,1)
hold on
h1 = plot(kcorrCompt_array, JK_array, '-k', "linewidth", 1)
##h1 = plot(kcorrCompt_array(ind_out_L), JK_array(ind_out_L), '-k', "linewidth", 1)
##h2 = plot(kcorrCompt_array(ind_out_R), JK_array(ind_out_R), '-k', "linewidth", 1)
h2 = plot(kcorrCompt_array(max_mu_Oskorbin_ind ), JK_array(max_mu_Oskorbin_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
hh = axis
hh = axis
xx = [kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(min(max_mu_Oskorbin_ind))]
yy = [hh(3) JK_array(min(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xx = [kcorrCompt_array(max(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))]
yy = [hh(3) JK_array(max(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xticks([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind)) ])
xticklabels([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))])
ylim([-0.28 -0.195])
ylabel('\it JK')
%
xx = [kcorrCompt_array(max_JK_array_ind) kcorrCompt_array(max_JK_array_ind)]
yy = [hh(3) JK_array(max_JK_array_ind)]
plot(xx, yy, '--r')
xticks([kcorrCompt_array(max_JK_array_ind)  ])
xticklabels([kcorrCompt_array(max_JK_array_ind) ])
ylim([-0.28 -0.195])
%
subplot(3,1,2)
hold on
h1 = plot(kcorrCompt_array, k_array, '-k', "linewidth", 1)
h2 = plot(kcorrCompt_array(max_mu_Oskorbin_ind ), k_array(max_mu_Oskorbin_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
xticklabels([])
yticklabels([])
yticks([1.7 1.8 1.9 2.0 ])
yticklabels([1.7 1.8 1.9 2.0 ])
ylim([1.7 2.1])
ylabel('\it K_{O}')
hh = axis
xx = [kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(min(max_mu_Oskorbin_ind))]
yy = [hh(3) k_array(min(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xx = [kcorrCompt_array(max(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))]
yy = [hh(3) k_array(max(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xticks([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind)) ])
xticklabels([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))])
%
subplot(3,1,3)
hold on
h1 = plot(kcorrCompt_array, max_mu_array, '-k', "linewidth", 1)
h2 = plot(kcorrCompt_array(max_mu_Oskorbin_ind ), max_mu_array(max_mu_Oskorbin_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
xticklabels([])
yticklabels([])
xticks([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind)) ])
xticklabels([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))])
yticks([17 18 ])
yticklabels([17 18 ])
ylim([16.6 18.5])
xlabel('\it k_{12}')
ylabel('\it max \mu')
set(gca, 'linewidth', 1);
set(gca, 'FontSize', 14);
hh = axis
xx = [kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(min(max_mu_Oskorbin_ind))]
yy = [hh(3) max_mu_array(min(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xx = [kcorrCompt_array(max(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))]
yy = [hh(3) max_mu_array(max(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
figure_name_out=strcat('PgammaPhTotalJKOskorbinMu', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

JK_arrayP = 0.5*(1+JK_array);
Kinv = 1/k_array;
max_mu_array01 =max_mu_array/numel(max_mu_array);

for jj = 1:numel(max_mu_array)
  T123 (jj)= JK_arrayP(jj)*Kinv(jj)*max_mu_array01(jj);
end

figure
hold on
plot(1:numel(max_mu_array), JK_arrayP)
plot(1:numel(max_mu_array), Kinv)
plot(1:numel(max_mu_array), max_mu_array01)
plot(1:numel(max_mu_array), 7*T123)

%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-10-10
% 2022-10-11
figure
subplot(4,1,1)
hold on
h1 = plot(kcorrCompt_array, JK_arrayP, '-k', "linewidth", 1)
##h1 = plot(kcorrCompt_array(ind_out_L), JK_array(ind_out_L), '-k', "linewidth", 1)
##h2 = plot(kcorrCompt_array(ind_out_R), JK_array(ind_out_R), '-k', "linewidth", 1)
h2 = plot(kcorrCompt_array(max_mu_Oskorbin_ind ), JK_arrayP(max_mu_Oskorbin_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
hh = axis
hh = axis
xx = [kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(min(max_mu_Oskorbin_ind))]
yy = [hh(3) JK_arrayP(min(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xx = [kcorrCompt_array(max(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))]
yy = [hh(3) JK_arrayP(max(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xticks([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind)) ])
xticklabels([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))])
ylim([0.36 0.4])
ylabel('\it JK')
%
xx = [kcorrCompt_array(max_JK_array_ind) kcorrCompt_array(max_JK_array_ind)]
yy = [hh(3) JK_arrayP(max_JK_array_ind)]
plot(xx, yy, '--r')
xticks([kcorrCompt_array(max_JK_array_ind)  ])
xticklabels([kcorrCompt_array(max_JK_array_ind) ])
%ylim([-0.28 -0.195])
%
subplot(4,1,2)
hold on
h1 = plot(kcorrCompt_array, Kinv, '-k', "linewidth", 1)
h2 = plot(kcorrCompt_array(max_mu_Oskorbin_ind ), Kinv(max_mu_Oskorbin_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
xticklabels([])
yticklabels([])
yticks([1.7 1.8 1.9 2.0 ])
yticklabels([1.7 1.8 1.9 2.0 ])
ylim([0.45 0.6])
ylabel('\it K_{O}')
hh = axis
xx = [kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(min(max_mu_Oskorbin_ind))]
yy = [hh(3) Kinv(min(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xx = [kcorrCompt_array(max(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))]
yy = [hh(3) Kinv(max(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xticks([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind)) ])
xticklabels([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))])
%
subplot(4,1,3)
hold on
h1 = plot(kcorrCompt_array, max_mu_array01, '-k', "linewidth", 1)
h2 = plot(kcorrCompt_array(max_mu_Oskorbin_ind ), max_mu_array01(max_mu_Oskorbin_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
xticklabels([])
yticklabels([])
xticks([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind)) ])
xticklabels([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))])
yticks([0.4 0.45 ])
yticklabels([0.4 0.45])
ylim([0.38 0.48])
%xlabel('\it k_{12}')
ylabel('\it max \mu')
set(gca, 'linewidth', 1);
set(gca, 'FontSize', 14);
hh = axis
xx = [kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(min(max_mu_Oskorbin_ind))]
yy = [hh(3) max_mu_array01(min(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xx = [kcorrCompt_array(max(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))]
yy = [hh(3) max_mu_array01(max(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')

subplot(4,1,4)
hold on
h1 = plot(kcorrCompt_array, T123, '-k', "linewidth", 1)
h2 = plot(kcorrCompt_array(max_mu_Oskorbin_ind ), T123(max_mu_Oskorbin_ind ), '-r', "linewidth", 1)
box('off')
set(h2, 'linewidth',1.5);
set(gca, 'FontSize', 14);
xticklabels([])
yticklabels([])
yticks([0.09 0.10 ])
yticklabels([0.09 0.10])
ylim([0.079 0.101])
%xlabel('\it k_{12}')
ylabel('\it T')
set(gca, 'linewidth', 1);
set(gca, 'FontSize', 14);
hh = axis
xx = [kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(min(max_mu_Oskorbin_ind))]
yy = [hh(3) T123(min(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xx = [kcorrCompt_array(max(max_mu_Oskorbin_ind)) kcorrCompt_array(max(max_mu_Oskorbin_ind))]
yy = [hh(3) T123(max(max_mu_Oskorbin_ind))]
plot(xx, yy, '--r')
xx = [kcorrCompt_array(max_JK_array_ind) kcorrCompt_array(max_JK_array_ind)]
yy = [hh(3) T123(max_JK_array_ind)]
plot(xx, yy, '--r')
xticks([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max_JK_array_ind) kcorrCompt_array(max(max_mu_Oskorbin_ind)) ])
xticklabels([kcorrCompt_array(min(max_mu_Oskorbin_ind)) kcorrCompt_array(max_JK_array_ind) kcorrCompt_array(max(max_mu_Oskorbin_ind))])

figure_name_out=strcat('PgammaPhTotalT', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%  PLOT 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
