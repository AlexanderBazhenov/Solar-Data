% 2022-10-04
clear all
close all

dirroot ='e:\Users\Public\Documents\ST\2022\T\'
dirData = 'e:\Users\Public\Documents\ST\2022\T\Solar\'
% 2022-04-14
dirroot = 'd:\Data\ST\2022\T\'
dirData = 'd:\Data\ST\2022\T\Solar\'


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get Data
cd(dirData), pwd
%
LIST1 = dir('Ch1*.csv');
LIST1.name;
LIST2 = dir('Ch2*.csv');
LIST2.name;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Octave interval
pkg load interval
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(dirroot), pwd
addpath(dirData)
addpath(dirroot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
FN1='Ch1_800nm_0.03.csv'
FN2='Ch2_800nm_0.03.csv'
%
FN1='Ch1_800nm_0.04.csv'
FN2='Ch2_800nm_0.04.csv'
%
FN1='Ch1_800nm_0.2.csv'
FN2='Ch2_800nm_0.2.csv'
%
FN1='Ch1_800nm_0.23mm.csv'
FN2='Ch2_800nm_0.23mm.csv'
%
FN1='Ch1_800nm_2mm.csv'
FN2='Ch2_800nm_2mm.csv'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-09-17
ii=13
% 2022-10-03
ii = 18
ii = 21 % 800nm_0.23mm
FN1=getfield(LIST1,{ii},'name');
FN2=getfield(LIST2,{ii},'name');
%
  epsilon0 = 10^(-4)
%  2022-12-27
%  epsilon0 = 2 *10^(-4)
  stepR = 0.01
[x1, x2, FNstr, Lambdastr, Threadstr] = getSolar2 (FN1, FN2);
input1 = x1(:,1);
input2 = x2(:,1);

% 2022-04-15
figure
hold on
h1=plot(x1(:,1),'.b')
h2=plot(x2(:,1),'.r')
lgd = legend (" Ch1 ", " Ch2 ");
set(lgd, 'location', "north");
set(lgd, 'location', "east");
text( 10, 0.99*max(x2(:,1)), FNstr, 'fontsize', 14 )
set(gca, 'fontsize', 14);
xlabel('n')
ylabel('mV')
figure_name_out=strcat(FNstr,'Raw', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% 2022-12-22
 figure;
 hold on
 nums = 1:numel(x2(:,1));
 [hax, h1, h2] = plotyy (nums, x1(:,1), nums, x2(:,1));
set (hax(1), 'fontsize', 14);
set (hax(1), "ycolor", [0 0 0]);
% set (hax, 'color', 'none');
set (hax(2), 'fontsize', 14);
set (hax(2), "ycolor", [0 0 0]);
xticks([50 100 150 200])
xticklabels([50 100 150 200])
xlabel('\it n')
set ([h1], "color", [1 0 0]);
set ([h2], "color", [0 0 1]);
lgd = legend (" Ch1 ", " Ch2 ");
set(lgd, 'location', "north");
set(lgd, 'fontsize', 14);
figure_name_out=strcat(FNstr,'Raw','2Axis', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Raw Data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% naive plot
R21=x2(:,1)./x1(:,1);
figure
hold on
h1=plot(R21,'.b')
set(gca, 'fontsize', 14);
text( 100, max(R21(:,1)), FNstr, 'fontsize', 14 )
xlabel('N')
ylabel('R21')
% 2022-09-17
ROIc = 50:150;
h2 = plot(ROIc, R21(ROIc),'.r')
figure_name_out=strcat(FNstr,'Ratio', '.png')
print('-dpng', '-r300', figure_name_out), pwd
maxR21 = max(R21(:,1))
minR21 = min(R21(:,1))
maxR21ROIc = max(R21(ROIc,1))
minR21ROIc = min(R21(ROIc,1))

% naive interlval plot
R21=x2(:,1)./x1(:,1);
figure
hold on
h1=plot(R21,'.b')
set(gca, 'fontsize', 14);
text( 100, max(R21(:,1)), FNstr, 'fontsize', 14 )
xlabel('N')
ylabel('R21')
% 2022-09-17
ROIc = 50:150;
h2 = plot(ROIc, R21(ROIc),'.r')
figure_name_out=strcat(FNstr,'Ratio', '.png')
print('-dpng', '-r300', figure_name_out), pwd
maxR21 = max(R21(:,1))
minR21 = min(R21(:,1))
maxR21ROIc = max(R21(ROIc,1))
minR21ROIc = min(R21(ROIc,1))

% naive Hist
figure
hist(R21(ROIc), 20)
hold on
set(gca, 'fontsize', 14);
xlabel('Rel. units')
ylabel('N')
text( 0.5*(maxR21ROIc + minR21ROIc), 1.05*max(hist(R21(ROIc), 20)), strcat('R21-ROIc', FNstr), 'fontsize', 14 )
dR = ( maxR21ROIc - minR21ROIc )/10
xlim([minR21ROIc-dR maxR21ROIc+dR])
%h2 = hist(R21(ROIc), 20)
%set(h2, 'color', 'r')
figure_name_out=strcat(FNstr,'HIST','ROIc', '.png')
print('-dpng', '-r300', figure_name_out), pwd

figure
hist(R21, 20)
hold on
set(gca, 'fontsize', 14);
xlabel('Rel. units')
ylabel('N')
text( 0.5*(maxR21 + minR21), 1.1*max(hist(R21, 20)), strcat('R21-', FNstr), 'fontsize', 14 )
dR = ( maxR21 - minR21 )/10
xlim([minR21-dR maxR21+dR])
%h2 = hist(R21(ROIc), 20)
%set(h2, 'color', 'r')
figure_name_out=strcat(FNstr,'HIST', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% /Vislual analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-12-23
input1 = x1(:,1);
xx1 = 1:length(input1);
xx = 1:length(input1)';
%epsilon0 = 10^(-4)
epsilon = epsilon0 * ones(length(input1),1);
SolarPair2errorbars
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% TRY

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%[tau, w, yint] = LinearModel (input1, epsilon0);
[tau1, w1, yint1] = DataLinearModelZ (input1, epsilon0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ChNo = 1
xx = 1:length(input1)';
rady1 = rad(yint1);
xx1 = 1:length(input1);
figure
h1=errorbar (xx1, mid(yint1), rady1, ".b");
set(h1,"linewidth",1)
set(h1,"markersize",10)
set(gca, 'fontsize', 14);
box('off')
%xp=x;
yp=tau1(1)+tau1(2)*xx;
hold on
h2=plot(xx1, yp, '-r')
set(gca, 'fontsize', 14)
%axis('tight')
haxis=axis
%xlim([haxis(1)-1 haxis(2)+1] )
ylim([haxis(3) haxis(4)] )
xlim([-5 length(mid(yint1))+5])
xlabel('n')
ylabel('mV')
%
text( 10, 0.99*max(yp), FNstr, 'fontsize', 14 )
ChText = strcat('Ch',num2str(ChNo))
text( 10, 0.95*max(yp), ChText, 'fontsize', 14 )
%
title_str=strcat('L1optimization', ' FN=',FN1)
title_str=strrep(title_str,' ','')
title_str=strrep(title_str,'_','-')
title_str=strrep(title_str,'.txt','')
title(title_str)
figure_name_out=strcat(title_str, '.png')
print('-dpng', '-r300', figure_name_out), pwd



figure
hist(w1)
xlabel('w')
ylabel('N')
box('off')
set(gca, 'fontsize', 14)
title_str=strcat('L1optimization-w', ' FN=',FN1)
title_str=strrep(title_str,' ','')
title_str=strrep(title_str,'_','-')
title_str=strrep(title_str,'.txt','')
title(title_str)
figure_name_out=strcat(title_str, '.png')
cd(dirData), pwd
print('-dpng', '-r300', figure_name_out), pwd

% remove linear drift
cd(dirroot)

yynew1 = yint1 -tau1(2)*xx1';
figure
h1=errorbar (xx, mid(yynew1), rady1, ".b");
set(h1,"linewidth",1)
set(h1,"markersize",10)
set(gca, 'fontsize', 14);
box('off')
set(gca, 'fontsize', 14)
haxis=axis
ylim([haxis(3) haxis(4)] )
xlim([1-0.5 length(yint1)+0.5])
%
title_str=strcat('Remove linear drift', ' FN=',FN1)
title_str=strrep(title_str,' ','')
title_str=strrep(title_str,'_','-')
title(title_str)
figure_name_out=strcat(title_str, '.png')
cd(dirData), pwd
print('-dpng', '-r300', figure_name_out), pwd
cd(dirroot)
% /remove linear drift

%%%%%%%%%%%%%%%%% NEW PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-10-03
ticklength = 0.5
lineWid = 1.5
Colors_special
colornow = OxfordBlue
figure
hold on
out = ErrorVert (xx, mid(yynew1), rady1, colornow, ticklength,  lineWid)
box('off')
set(gca, 'fontsize', 14)
haxis=axis
ylim([haxis(3) haxis(4)] )
xlim([1-0.5 length(yint1)+0.5])
pbaspect([2 1 1])
figure_name_out=strcat(FNstr,'Corrected','Ch1','Oskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% Equal Errors
radyy = epsilon0*ones(length(xx),1);
figure
hold on
out = ErrorVert (xx, mid(yynew1), radyy, colornow, ticklength,  lineWid)
box('off')
set(gca, 'fontsize', 14)
haxis=axis
ylim([haxis(3) haxis(4)] )
xlim([1-0.5 length(yint1)+0.5])
pbaspect([2 1 1])
figure_name_out=strcat(FNstr,'Corrected','Ch1','Eq', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%% /NEW PLOT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
hist(mid(yynew1), 20)
set(gca, 'fontsize', 14);
xlabel('Rel. units')
ylabel('N')
title_str=strcat('Remove linear drift', ' HIST', ' FN=',FN1)
title_str=strrep(title_str,' ','')
title_str=strrep(title_str,'_','-')
title(title_str)
figure_name_out=strcat(FN1, ' RemLin',' HIST', '.png')
cd(dirData), pwd
print('-dpng', '-r300', figure_name_out), pwd
cd(dirprocessing)
% /remove linear drift


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% REF

yy = x2(:,1);
epsilon = epsilon0 * ones(length(yy),1);
inf_b = yy - epsilon;
sup_b = yy + epsilon;
[tau2, w2, yint2] = DataLinearModelZ (input2, epsilon0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear programming
%[tau, w] = L_1_minimization(X, inf_b, sup_b);
%
yynew2 =  yint2 -tau2(2)*xx';
rady2 = rad(yynew2);
figure
h1=errorbar (xx, mid(yynew2), rady2, ".b");

figure
hist(mid(yynew2), 20)
set(gca, 'fontsize', 14);
xlabel('Rel. units')
ylabel('N')
title_str=strcat('Remove linear drift', ' HIST', ' FN=',FN2str)
title(title_str)
figure_name_out=strcat(FN2str, ' RemLin',' HIST', '.png')
cd(dirData), pwd
print('-dpng', '-r300', figure_name_out), pwd
cd(dirprocessing)
cd(dirroot)

% 2022-10-04
ChNo = 2
rady2 = rad(yint2);
figure
hold on
out = ErrorVert (xx, yynew2, rady2, colornow, ticklength,  lineWid)
box('off')
set(gca, 'fontsize', 14)
haxis=axis
ylim([haxis(3) haxis(4)] )
xlim([1-0.5 length(yint1)+0.5])
pbaspect([2 1 1])
figure_name_out=strcat(FNstr,'Corrected','Ch2','Oskorbin', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% Equal Errors
radyy = epsilon0*ones(length(xx),1);
figure
hold on
out = ErrorVert (xx, yynew2, radyy, colornow, ticklength,  lineWid)
box('off')
set(gca, 'fontsize', 14)
haxis=axis
ylim([haxis(3) haxis(4)] )
xlim([1-0.5 length(yint1)+0.5])
pbaspect([2 1 1])
figure_name_out=strcat(FNstr,'Corrected','Ch2','Eq', '.png')
print('-dpng', '-r300', figure_name_out), pwd

% 2022-12-21
figure
hold on
h1=errorbar (xx, mid(yynew1)+0.03, rady1, ".b");
h2=errorbar (xx, mid(yynew2), rady2, ".r");
set(h1,"linewidth",1)
set(h1,"markersize",10)
set(gca, 'fontsize', 14);
box('off')
set(gca, 'fontsize', 14)
haxis=axis
ylim([haxis(3) haxis(4)] )
xlim([1-0.5 length(yint1)+0.5])
yticklabels([])
figure_name_out=strcat(FNstr,'Corrected','Ch1Ch2','Eq', '.png')
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-12-22
[tau1, w1, yint1] = DataLinearModelZ (input1, epsilon0);
yynew1 = yint1 -tau1(2)*xx1';
S = yynew1;
[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S), sup(S)]);
outer1 = infsup( min(inf(S)), max(sup(S)) )
inner1 = infsup( max(inf(S)), min(sup(S)) )
wid(outer1)
%
[tau2, w2, yint2] = DataLinearModelZ (input2, epsilon0);
yynew2 = yint2 -tau2(2)*xx1';
S = yynew2;
[mode2, modefreq2, freqs2, Ss2] = imodeR([inf(S), sup(S)]);
outer2 = infsup( min(inf(S)), max(sup(S)) )
inner2 = infsup( max(inf(S)), min(sup(S)) )
wid(outer2)
R21outer =  outer2 / outer1
wid(R21outer)
R21inner =  inner2 / inner1

inner2/inner1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

rady1 = rad(yint1);
rady2 = rad(yint2);
% 2022-12-23
SolarPair2Consterrorbars

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R21min=min(yynew2)/max(yynew1)
R21max=max(yynew2)/min(yynew1)
% Find Jaccard array
cd(dirroot)
JK_array = []
k_array = []
m_array = []
R21est=1.065:0.00001:1.066;
%
for R21now=R21est
  X = [ yynew2; R21now*yynew1];
 % Xint=midrad(X, epsilon0  * ones(length(X),1));
  JK_now =jaccardKRSet(X);
  JK_array = [JK_array JK_now];
  [oskorbin_center_k, k_now] = estimate_uncertainty_center(X);
  k_array = [k_array k_now];
  [mode_now, modefreq_now, freqs_now, Ss_now] = imodeR([inf(X), sup(X)]);
  m_array = [m_array modefreq_now];
end
[max_JK_array max_JK_array_ind] = max(JK_array)
R21opt = R21est(max_JK_array_ind)

% Plot JK Oskorbin Mu
PlotJKOskorbinMuSolar
figure_name_out=strcat(FNstr,'JKOskorbinMu', '.png')
print('-dpng', '-r300', figure_name_out), pwd

JK_arrayP = (1 + JK_array)/2;
Kinv = 1./k_array;
max_mu_array01 = m_array/length(X);
for jj=1:length(R21est)
  T123 (jj) = JK_arrayP(jj)*Kinv(jj)*max_mu_array01(jj);
end
% Plot Tot = JK * 1/Oskorbin * Mu
PlotJKOskorbinMuTotSolar
figure_name_out=strcat(FNstr,'SolarTotalT', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Xest = [ yynew2; R21opt*yynew1];
S = Xest;
[mode, modefreq, freqs, Ss] = imodeR([inf(S), sup(S)]);

S = yynew2;
[mode2, modefreq2, freqs2, Ss2] = imodeR([inf(S), sup(S)]);

S = R21opt*yynew1;
[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S), sup(S)]);


figure
plot(Ss(1:end-1), freqs)
hold on
plot(Ss2(1:end-1), freqs2)
plot(Ss1(1:end-1), freqs1)
legend('Total', 'X2', 'R*X1')
set(gca, 'fontsize', 14)
ylim( [0 400])
xlim([min(Ss)-epsilon0/4 max(Ss)+epsilon0/4 ] )
xlabel('\it R')
xlabel('')
ylabel('\it \mu')
pbaspect([2 1 1])
figure_name_out=strcat(FNstr,'Freq','Combined', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all
Xest0 = [ yynew2; yynew1];
Xest2 = [ yynew2; (R21opt-0.00435)*yynew1];
jaccardKRSet(Xest2)
figure
hold on
h0=errorbar (1:length(Xest0), mid(Xest0), rad(Xest0),".k");
h1=errorbar (1:length(Xest), mid(Xest), rad(Xest),".b");
h2=errorbar (1:length(Xest2), mid(Xest2), rad(Xest2),".r");
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure
h = plot(R21est, JK_array)
set(gca, 'fontsize', 14);
axis('tight')
box('off')
title_str=strcat('Jaccard vs 21Ratio',  ' FN=',FNstr)
title(title_str)
figure_name_out=strcat(FNstr, ' RemLin',' HIST', '.png')
cd(dirData), pwd
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-08-25
% generic figure
figure
h = plot(R21est, JK_array)
set(gca, 'fontsize', 14);
axis('tight')
box('off')
set(gca, 'fontsize', 14);

set(gca, 'linewidth', 1);
set(h, 'linewidth', 2);
hold on
xx = [min(R21est) max(R21est)]
yy = [-1.2*abs(max(JK_array )) -1.2*abs(max(JK_array)) ]
plot(xx, yy, '--k')
indmore = find(JK_array > -1.2*abs(max(JK_array )) )
xx = [ min(indmore) max(indmore) ]
yy = [ JK_array(min(indmore)) JK_array(max(indmore)) ]
plot(xx, yy, 'or')
xx = [ R21est(min(indmore)) R21est(min(indmore)) ]
yy = [ min(JK_array) JK_array(min(indmore)) ]
plot(xx, yy, '--k')
xx = [ R21est(max(indmore)) R21est(max(indmore)) ]
yy = [ min(JK_array) JK_array(max(indmore)) ]
plot(xx, yy, '--k')
[ JKmax, indmax ] = max(JK_array)
xx = [ min(R21est) R21est(indmax) ]
yy = [ JKmax JKmax ]
plot(xx, yy, '--k')
ylim([ min(JK_array) 0.95*JKmax] )

xticks([])
yticks([])
title('')
xx = [ R21est(min(indmore)) R21est(max(indmore)) ]
yy = [ min(JK_array) min(JK_array) ]
h1 = plot(xx, yy, '-r')
set(h1, 'linewidth', 4);

xx = [ min(R21est) min(R21est)  ]
yy = [ -1.2*abs(max(JK_array )) max(JK_array) ]
h2 = plot(xx, yy, '-r')
set(h2, 'linewidth', 4);

figure_name_out=strcat('JKvsParam', '.png')
cd(dirroot), pwd
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cd(dirroot), pwd
  X = [ yynew2; maxR21*yynew1];
  XY = midrad(X, epsilon0  * ones(length(X),1));
 [mode, mu_array, max_mu, mode_ind, c_array, C, multi]= modeIR4(XY)
[oskorbin_center, k] = estimate_uncertainty_center(XY)
figure
h = plot(mu_array)
set(h, 'linewidth', 2)
box('off')
ylabel('\it \mu')
xlabel('c intervals')
set(gca, 'fontsize', 14);
title_str=strcat(' mu array -',  FNstr)
title(title_str)
figure_name_out=strcat(FNstr, ' mu array', '.png')
cd(dirData), pwd
print('-dpng', '-r300', figure_name_out), pwd

 [mode1, mu_array1, max_mu1, mode_ind1, c_array1, C1, multi1]= modeIR4(yynew1)
 [oskorbin_center1, k1] = estimate_uncertainty_center(yynew1)
 figure
h = plot(mu_array1)
set(h, 'linewidth', 2)
box('off')
ylabel('\it \mu')
xlabel('c intervals')
set(gca, 'fontsize', 14);
title_str=strcat(' mu array Ch1-',  FNstr)
title(title_str)
figure_name_out=strcat(FNstr, ' mu array Ch1', '.png')
cd(dirData), pwd
print('-dpng', '-r300', figure_name_out), pwd

