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



# Значения независимой переменной x и
# соответствующие им интервальные наблюдения откликовой переменной y
##x = [   3;         6;        9;       12;       15;       18  ];
##y = [ 12 17;     10 13;    13 18;    14 18;    19 24;   16 19 ];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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

epsilon0 = 10^(-4)
epsilon = epsilon0 * ones(length(input1),1);

%x = xx1';
for ii = 1:length(input1)
  yy(ii,1) = input1(ii)-epsilon0;
  yy(ii,2) = input1(ii)+epsilon0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xx = xx1';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = xx;
y = yy;

# График исходных данных
figure
ir_scatterlu(xx, y);
title('Исходные данные')
xlabel('x');
ylabel('y');

# 1. C заданными точками излома
p = [ x(100)];  # точки излома (могут быть отличными от точек выборки)
p = [  x(50); x(100); x(150)];
ROIp=[ 50 150 ]
ROIp=[ 20 180 ]
ROIp=[ 25 175 ]
ROIp=[ 50 100 150 ]
ROIp=[ 25 50 75 100 125 150 175 ]
ROIp=[ 12 25 37 50 62 75 87 100 112 125 137 150 162 175 ]

##p = x(ROIp)
##SolarPiecewiceManualDo

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ROIp=[ 50 100 150 ]
ROIp=[ 12 25 37 50 62 75 87 100 112 125 137 150 162 175 ]
clear y
clear x

n = numel(ROIp)
y(:,1) = yy(ROIp,1)
y(:,2) = yy(ROIp,2)
x = xx(ROIp)
x = 1:numel(ROIp)
x = x'
X = [ ones(n,1) (1:n)']
figure
ir_scatterlu(X, y);
title('Исходные данные')
xlabel('x');
ylabel('y');

# 2. C автоматически выбираемыми точками излома
C = numel(ROIp)-2;          #  максимально допустимое количество выбираемых точек излома
C = 2
delta = 10^(-4) #  порог значимости коэффициентов кусочно-линейной части
##delta = delta / 2; #  порог значимости коэффициентов кусочно-линейной части
##delta = delta * 1.1; #  порог значимости коэффициентов кусочно-линейной части
[beta, exitcode, p] = ir_estimatepwa(x,y,C,delta);
exitcode

if exitcode == 0
    disp('Коэффициенты регрессии:')
    disp(beta)
    disp('Точки излома:')
    disp(p)

    disp('Коэффициенты прямых-кусков:')
    gamma = ir_pieces(beta, p)
else
    disp('Ошибка вычисления коэффициентов!')
    return
end

# Вычисление интервальных прогнозных значений кусочно-линейной регрессии
# в точках x
[yp] = ir_predictpw(x, beta, p);
##disp('Прогнозные значения регрессии:')
##    disp(yp);

# График исходных значений, прогнозных значений и
# множества допустимых моделей
figure
hold on
ir_plotmodelsetpw(x, beta, p);  # График множества допустимых моделей
ir_scatterlu(x, y,  'k.');      # Исходные данные
%ir_scatterlu(x, yp, 'm.');      # Предсказания в исходных точках x
title({'Множество допустимых кусочно-линейных моделей';...
       ['с изломами в точках: ' sprintf('%1.1f ',p)]});
xlabel('x');
ylabel('y');


Pstr = num2str(p')
PstrNB = strrep(Pstr,' ','')
figure_name_out=strcat(FNstr,'Piecewice',PstrNB, '.png')
print('-dpng', '-r300', figure_name_out), pwd

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 2023-03-21
Angle_points = [ 25 176 ]
%
ROI_array = zeros(numel(Angle_points)+1,2)
ROI_array (1,1) = 1
ROI_array (1,2) = Angle_points(1)-1
ROI_array (2,1) = Angle_points(1)
ROI_array (2,2) = Angle_points(2)-1
ROI_array (end,1) =  Angle_points(2)
ROI_array (end,2) = length(input1)


figure
hold on
for ii = 1:length(ROI_array)
 % clear yp;
  ROI_now = ROI_array(ii,1):  ROI_array(ii,2);
  data_now = input1(ROI_now);
  in_now = ROI_now;
 % [tau1, w1, yint] = DataLinearModelZ (data_now, epsilon0);
  h1=errorbar (in_now, data_now, epsilon(ROI_now), ".k");
%
%  yp= tau1(1) +tau1(2)*(in_now-in_now(1)+1);


h2=plot(in_now, yp(in_now,1), '-r')
h3=plot(in_now, yp(in_now,2), '-r')

  end
set(gca, 'fontsize', 14)
ylim([min(input1)-epsilon0 max(input1)+epsilon0 ] )
xlim([1 200])

% Remove  Multi Linear
ypmid = (yp(:,1)+yp(:,2))/2;
yynew=input1-ypmid+input1(1);
y_now=yynew;
figure
h1=errorbar (x, y_now, epsilon, ".b");
ylim([min(y_now)-epsilon0 max(y_now)+epsilon0 ] )
xlim([1 200])
% JK
S1 = midrad(y_now, epsilon );
JK_now =jaccardKRSet(S1);
[mode1, modefreq1, freqs1, Ss1] = imodeR([inf(S1), sup(S1)]);