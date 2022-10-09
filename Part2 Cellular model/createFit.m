function [fitresult, gof] = createFit(x, y)
%CREATEFIT(X,Y)
%  创建一个拟合。
%
%  要进行 'The fitting of experimental data' 拟合的数据:
%      X 输入: x
%      Y 输出: y
%  输出:
%      fitresult: 表示拟合的拟合对象。
%      gof: 带有拟合优度信息的结构体。
%
%  另请参阅 FIT, CFIT, SFIT.

%  由 MATLAB 于 09-Oct-2022 14:15:35 自动生成


%% 拟合: 'The fitting of experimental data'。
[xData, yData] = prepareCurveData( x, y );

% 设置 fittype 和选项。
ft = fittype( '(n*(x^(n))*17.5618)/(k+x^(n))', 'independent', 'x', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.StartPoint = [0.678735154857773 0.757740130578333];

% 对数据进行模型拟合。
[fitresult, gof] = fit( xData, yData, ft, opts );

% 绘制数据拟合图。
figure( 'Name', 'The fitting of experimental data');
h = plot( fitresult, xData, yData,"*");
legend( h, 'Experimental data', 'The fitting of experimental data', 'Location', 'Best', 'Interpreter', 'none' );
% 为坐标区加标签
title('The fitting of experimental data','fontsize', 16);
set(gca, 'linewidth', 1.1, 'fontsize', 16, 'fontname', 'times');
xlabel( 'The concentration of unbound silver ions (μM)', 'Interpreter', 'none',"FontSize",16 );
ylabel( 'The concentration of absorbed AgNP(silver ions) (μM)', 'Interpreter', 'none',"FontSize",16 );
grid on


