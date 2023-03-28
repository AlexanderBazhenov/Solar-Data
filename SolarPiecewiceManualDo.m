% 2023-02-23

[beta, exitcode] = ir_estimatepw(x,y,p);
exitcode

# 2. C автоматически выбираемыми точками излома
##C = 3;          #  максимально допустимое количество выбираемых точек излома
##C = 3;          #  максимально допустимое количество выбираемых точек излома
##delta = 0.0001; #  порог значимости коэффициентов кусочно-линейной части
##delta = delta / 2; #  порог значимости коэффициентов кусочно-линейной части
##delta = delta * 1.1; #  порог значимости коэффициентов кусочно-линейной части
##[beta, exitcode, p] = ir_estimatepwa(x,y,C,delta);
##exitcode

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
%title({'Множество допустимых кусочно-линейных моделей';...
%       ['с изломами в точках: ' sprintf('%1.1f ',p)]});
title_str=strcat('Angle points: ', sprintf('%2.0f ',p));
title(title_str);
xlabel('x');
ylabel('y');


Pstr = num2str(p')
PstrNB = strrep(Pstr,' ','')
figure_name_out=strcat(FNstr,'Piecewice',PstrNB, '.png')
print('-dpng', '-r300', figure_name_out), pwd

