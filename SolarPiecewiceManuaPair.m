% 2023-03-21
% 2023-02-23
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

ii = 21 % 800nm_0.23mm
FN1=getfield(LIST1,{ii},'name');
FN2=getfield(LIST2,{ii},'name');
%  stepR = 0.01
[x1, x2, FNstr, Lambdastr, Threadstr] = getSolar2 (FN1, FN2);
input1 = x1(:,1);
input2 = x2(:,1);
xx1 = 1:length(input1);

epsilon0 = 10^(-4)
epsilon = epsilon0 * ones(length(input1),1);

ROIp=[ 25 175 ]
ROIp=[ 20 180 ]
xx = xx1';

%x = xx1';
for ii = 1:length(input1)
  yy(ii,1) = input1(ii)-epsilon0;
  yy(ii,2) = input1(ii)+epsilon0;
end
x = xx;
y = yy;
p = x(ROIp)

##SolarPiecewiceManualDo
[beta, exitcode] = ir_estimatepw(x,y,p);
[yp] = ir_predictpw(x, beta, p);
% Remove  Multi Linear
ypmid1 = (yp(:,1)+yp(:,2))/2;
yynew1=input1-ypmid1+input1(1);
y_now=yynew1;
% 
S1 = midrad(y_now, epsilon );
%
for ii = 1:length(input1)
  yy(ii,1) = input2(ii)-epsilon0;
  yy(ii,2) = input2(ii)+epsilon0;
end
y = yy;
##SolarPiecewiceManualDo
[beta, exitcode] = ir_estimatepw(x,y,p);
[yp] = ir_predictpw(x, beta, p);
% Remove  Multi Linear
ypmid2 = (yp(:,1)+yp(:,2))/2;
yynew2=input2-ypmid2+input2(1);
y_now=yynew2;
S2 = midrad(y_now, epsilon );

%%%%%%%%%%%%%%%%%%%%%% LINEAR MODEL TREND %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S1), sup(S1)]);
outer1 = infsup( min(inf(S1)), max(sup(S1)) )
inner1 = infsup( max(inf(S1)), min(sup(S1)) )
wid(inner1)
wid(outer1)
%
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
SolarPair2Consterrorbars
figure_name_out=strcat(FNstr,'InteConst','Eps=',num2str(epsilon0), '.png')
print('-dpng', '-r300', figure_name_out), pwd
SolarPair2ModeS1S2
figure_name_out=strcat(FNstr,'Freq','S1S2w0', '.png')
print('-dpng', '-r300', figure_name_out), pwd

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
%  [oskorbin_center_k, k_now] = estimate_uncertainty_center(X);
    [oskorbin_center_k, k_now] = estimate_uncertainty_center1(X);
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


