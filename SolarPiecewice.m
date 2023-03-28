% 2023-02-22

# Значения независимой переменной x и
# соответствующие им интервальные наблюдения откликовой переменной y
##x = [   3;         6;        9;       12;       15;       18  ];
##y = [ 12 17;     10 13;    13 18;    14 18;    19 24;   16 19 ];


input1 = x1(:,1);
input2 = x2(:,1);
xx1 = 1:length(input1);
epsilon0 = 10^(-4)
epsilon = epsilon0 * ones(length(input1),1);

x = xx1';
for ii = 1:length(input1)
  y(ii,1) = input1(ii)-epsilon0;
  y(ii,2) = input1(ii)+epsilon0;
end
# График исходных данных
figure
ir_scatterlu(x, y);
title('Исходные данные')
xlabel('x');
ylabel('y');

# 1. C заданными точками излома
p = [ x(50); x(150) ];  # точки излома (могут быть отличными от точек выборки)
[beta, exitcode] = ir_estimatepw(x,y,p);
p = [ x(25); x(175) ];  # точки излома (могут быть отличными от точек выборки)
[beta, exitcode] = ir_estimatepw(x,y,p);

# 2. C автоматически выбираемыми точками излома
C = 3;          #  максимально допустимое количество выбираемых точек излома
C = 2;          #  максимально допустимое количество выбираемых точек излома
delta = 0.0001; #  порог значимости коэффициентов кусочно-линейной части
delta = delta / 2; #  порог значимости коэффициентов кусочно-линейной части
delta = delta * 1.1; #  порог значимости коэффициентов кусочно-линейной части
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
disp('Прогнозные значения регрессии:')
    disp(yp)

# График исходных значений, прогнозных значений и
# множества допустимых моделей
figure
hold on
ir_plotmodelsetpw(x, beta, p);  # График множества допустимых моделей
ir_scatterlu(x, y,  'k.');      # Исходные данные
ir_scatterlu(x, yp, 'm.');      # Предсказания в исходных точках x
title({'Множество допустимых кусочно-линейных моделей';...
       ['с изломами в точках: ' sprintf('%1.1f ',p)]});
xlabel('x');
ylabel('y');


Pstr = num2str(p')
PstrNB = strrep(Pstr,' ','')
figure_name_out=strcat(FNstr,'Piecewice',PstrNB, '.png')

figure_name_out=strcat(FNstr,'Piecewice','Auto', '.png')
print('-dpng', '-r300', figure_name_out), pwd

