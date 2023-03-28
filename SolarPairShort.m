% 2022-10-04
clear all
close all

dirroot ='e:\Users\Public\Documents\ST\2023\T\', dirData = 'e:\Users\Public\Documents\ST\2022\T\Solar\'
dirOld =  'e:\Users\Public\Documents\ST\2022\T\'
% 2022-04-14
dirroot = 'd:\Data\ST\2023\T\', dirData = 'd:\Data\ST\2022\T\Solar\'
dirOld =  'd:\Data\ST\2022\T\'
%
dirpiecewise = strcat(dirroot,'ir_piecewise')

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
addpath(dirOld)
addpath(dirpiecewise)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
% epsilon0 = 10 *10^(-4)
  stepR = 0.01
[x1, x2, FNstr, Lambdastr, Threadstr] = getSolar2 (FN1, FN2);
input1 = x1(:,1);
input2 = x2(:,1);
xx1 = 1:length(input1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2023-02-15
figure
h1=plot (xx1, input1, ".b");
hold on
h2=plot(xx1, input2,  ".r");
set(gca, 'fontsize', 14)
xlim([0 length(xx1)+1])
xlabel('\it n')
ylabel('\it mV')
pbaspect([1 1 1])
figure_name_out=strcat(FNstr,'Raw', '.png')
print('-dpng', '-r300', figure_name_out), pwd





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-12-23
input1 = x1(:,1);
xx1 = 1:length(input1);
% epsilon0 = 10^(-4)
epsilon = epsilon0 * ones(length(input1),1);
%%%%%%%%%%%%%%%%%% PLOT RAW DATA %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SolarPair2errorbars
figure_name_out=strcat(FNstr,'Inte','2Axis', '.png')
print('-dpng', '-r300', figure_name_out), pwd

X1 = midrad(input1, epsilon);
JK1raw = jaccardKRSet(X1)
X2 = midrad(input2, epsilon);
JK2raw = jaccardKRSet(X2)

%%%%%%%%%%%%%%%%%%%% PLOT Linear Models %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SolarPair2LinModels
figure_name_out=strcat(FNstr,'Inte','2Axis','LinMode', '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%% LINEAR MODEL TREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-12-22
% w > 1
[tau1, w1, yint1] = DataLinearModel (input1, epsilon0);
yynew1 = input1 -tau1(2)*xx1';
S1 = midrad(yynew1, epsilon );
[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S1), sup(S1)]);
[tau2, w2, yint2] = DataLinearModel(input2, epsilon0);
yynew2 = input2-tau2(2)*xx1';
S2 = midrad(yynew2, epsilon );
[mode2, modefreq2, freqs2, Ss2] = imodeR([inf(S2), sup(S2)]);
%
SolarPair2ModeS1S2
figure_name_out=strcat(FNstr,'Freq','S1S2w1', '.png')
print('-dpng', '-r300', figure_name_out), pwd
%
X1line = midrad(yynew1, epsilon);
JK1line = jaccardKRSet(X1line)
X2line = midrad(yynew2, epsilon);
JK2line= jaccardKRSet(X2line)


%%%%%%%%%%%%%%%%%%%%%% LINEAR MODEL TREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-12-27
% w > 0
[tau1, w1, yint] = DataLinearModelZ (input1, epsilon0);
%yynew1 = yint1 -tau1(2)*xx1';
% 2022-12-28
yynew1 = input1 -tau1(2)*xx1';
S1 = midrad(yynew1, epsilon );
%
[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S1), sup(S1)]);
outer1 = infsup( min(inf(S1)), max(sup(S1)) )
inner1 = infsup( max(inf(S1)), min(sup(S1)) )
wid(inner1)
wid(outer1)
%
% [tau2, w2, yint2] = DataLinearModel(input2, epsilon0);
[tau2, w2, yint2] = DataLinearModelZ(input2, epsilon0);
yynew2 = input2-tau2(2)*xx1';
S2 = midrad(yynew2, epsilon );
[mode2, modefreq2, freqs2, Ss2] = imodeR([inf(S2), sup(S2)]);
outer2 = infsup( min(inf(S2)), max(sup(S2)) )
inner2 = infsup( max(inf(S2)), min(sup(S2)) )
wid(inner2)
wid(outer2)
%
R21outer =  outer2 / outer1
wid(R21outer)
R21inner =  inner2 / inner1
wid(R21inner)

X1line = midrad(yynew1, epsilon);
JK1line = jaccardKRSet(X1line)
X2line = midrad(yynew2, epsilon);
JK2line= jaccardKRSet(X2line)

%%%%%%%%%%%%%%%%%%%%% LINEAR MODEL TREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2022-12-23
% 2022-12-28
SolarPair2Consterrorbars
% w>1
% figure_name_out=strcat(FNstr,'InteConst','w1','Eps=',num2str(epsilon0), '.png')
% w> 0
figure_name_out=strcat(FNstr,'InteConst','Eps=',num2str(epsilon0), '.png')
print('-dpng', '-r300', figure_name_out), pwd

SolarPair2w

%%%%%%%%%%%%%%%%%%%%%%%%%% Mode S1 S2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
##S = S2;
##[mode2, modefreq2, freqs2, Ss2] = imodeR([inf(S), sup(S)]);
##S = S1;
##[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S), sup(S)]);

SolarPair2ModeS1S2
figure_name_out=strcat(FNstr,'Freq','S1S2w0', '.png')
print('-dpng', '-r300', figure_name_out), pwd


%%%%%%%%%%%%%%%%% Multiliear Model ROI_array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 2022-01-07
%
Angle_points = [ 50 151 ]
% 2022-03-21
Angle_points = [ 25 176 ]
%
ROI_array = zeros(numel(Angle_points)+1,2)
ROI_array (1,1) = 1
ROI_array (1,2) = Angle_points(1)-1
ROI_array (2,1) = Angle_points(1)
ROI_array (2,2) = Angle_points(2)-1
ROI_array (end,1) =  Angle_points(2)
ROI_array (end,2) = length(input1)
%%%%%%%%%%%%%%%%%%%%%%%%%%%% ROI_array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Multi Linear
SolarPair2LinMultiModels
figure_name_out=strcat(FNstr,'MultiModel', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% Remove  Multi Linear
SolarPair2MultiModels
figure_name_out=strcat(FNstr,'MultiModelCorrected', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% Parameters
S1 = midrad(yynew1, epsilon );
S2 = midrad(yynew2, epsilon );
% Mode
S = S2;
[mode2, modefreq2, freqs2, Ss2] = imodeR([inf(S), sup(S)]);
S = S1;
[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S), sup(S)]);
SolarPair2ModeS1S2
figure_name_out=strcat(FNstr,'Freq','MultiModel', '.png')
print('-dpng', '-r300', figure_name_out), pwd
% Jaccard
JKMM1 = jaccardKRSet(S1)
JKMM2 = jaccardKRSet(S2)
%
outer1 = infsup( min(inf(S1)), max(sup(S1)) )
inner1 = infsup( max(inf(S1)), min(sup(S1)) )
wid(inner1)
wid(outer1)
outer2 = infsup( min(inf(S2)), max(sup(S2)) )
inner2 = infsup( max(inf(S2)), min(sup(S2)) )
wid(inner2)
wid(outer2)
%
R21outer =  outer2 / outer1
wid(R21outer)
R21inner =  inner2 / inner1
wid(R21inner)
%%%%%%%%%%%%%%%%%%%%% QUALITY FUNCTIONALS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R21min=min(yynew2)/max(yynew1)
R21max=max(yynew2)/min(yynew1)
% Find Jaccard array
cd(dirroot)
JK_array = []
k_array = []
m_array = []
R21est=1.065:0.00001:1.0665;
%
for R21now=R21est
  X = [ S2; R21now*S1];
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
figure_name_out=strcat(FNstr,'JKOskorbinMu','Eps=',num2str(epsilon0), 'AlphaLev=',num2str(alpha_lev) ,'.png')
print('-dpng', '-r300', figure_name_out), pwd

JK_arrayP = (1 + JK_array)/2;
Kinv = 1./k_array;
%Kinv = Kinv/max (Kinv);
max_mu_array01 = m_array/length(X);
for jj=1:length(R21est)
  T123 (jj) = JK_arrayP(jj)*Kinv(jj)*max_mu_array01(jj);
end
% Plot Tot = JK * 1/Oskorbin * Mu
[max_T123 max_T123_ind] = max(T123)
max_T123_ind_more = find(T123 > 0.99*max_T123)
% 2023-01-19
alpha_lev = 0.99 %0.99%
max_T123_ind_more = find (T123> alpha_lev*max_T123 )

Rint = infsup(R21est(min(max_T123_ind_more)), R21est(max(max_T123_ind_more)))

alpha_lev_mu = 0.985 %0.999%
PlotJKOskorbinMuTotSolar
R21est(min(alpha_set_mu))
R21est(max(alpha_set_mu))
Rout = infsup(R21est(min(alpha_set_mu)), R21est(max(alpha_set_mu)))
 figure_name_out=strcat(FNstr,'SolarTotalT','Eps=',num2str(epsilon0), 'AlphaLev=',num2str(alpha_lev) ,'.png')
print('-dpng', '-r300', figure_name_out), pwd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


jaccardKR(Rint, Rout )
