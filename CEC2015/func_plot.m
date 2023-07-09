% function function_plot
clear;close all
global initial_flag

% Xmin=[-100,-100,-100,-100,-100,-100,0,-32,-5,-5,-0.5,-pi,-3,-100,-5,-5,-5,-5,-5,-5,-5,-5,-5,-5,2];
% Xmax=[100,100,100,100,100,100,600,32,5,5,0.5,pi,1,100,5,5,5,5,5,5,5,5,5,5,5];

for func_num=5
    if func_num==1 x=-100:5:100;y=x; %[-100,100]
    elseif func_num==2 x=-100:2:100; y=x;%[-10,10]
    elseif func_num==3 x=-110:1:100; y=-100:1:110;%[-100,100]
    elseif func_num==4 x=-100:1:100; y=x;%[-100,100]
    elseif func_num==5 x=-100:1:100;y=x; %[-5,5]
    elseif func_num==6 x=78:0.05:82; y=-51.3:0.05:-47.3;%[-3,1]
    elseif func_num==7 x=-350:2:-250; y=-100:2:0;;%[-600,600]
    elseif func_num==8 x=-32:1:32; y=x;%[-32,32]
    elseif func_num==9 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==10 x=-5:0.1:5; y=x;%[-5,5]
    elseif func_num==11 x=-100:1:100; y=x;%[-0.5,0.5]
    elseif func_num==12 x=-3:0.1:3; y=x;%[-pi,pi]
    elseif func_num==13 x=-2:0.02:-1; y=x;%[-3,1]
    elseif func_num==14 x=-90:0.2:-50; y=-40:0.2:0;%[-100,100]
    elseif func_num==15 x=-5:0.1:5; y=x;%[-5,5]
    end

    initial_flag=0;
    L=length(x);
    L1=length(y);
    f=[];

    for i=1:L
        for j=1:L
            f(i,j)=benchmark_func([x(i),y(j)],func_num);
        end
    end


    figure(func_num);
    surfc(y,x,f);
%     plot(x,f);
%     contourf(y,x,f);%等高线图
%     colormap(bone);
    meshc(x,y,f);%三维曲面（浅色）a+等高线
%     surfl(x,y,f);
%     shading interp;
%     colormap(gray);
    hold on
    
    result = xlsread("D:/智能优化/DE变体相关文章与源码/DE_Framework_CEC2015_30D/收敛性分析/f4.xlsx");
    for i = 1:5:size(result,1)
        X = result(i,:);
%         if (-80 < X(1) && X(1) < -20) && (-80 < X(2) && X(2) < -20)
        
        plot3(X(2),X(1),benchmark_func([X(1),X(2)],func_num),'r*','linewidth',1.5);
%         end
        hold on
    end
    xlabel("X");
    ylabel("Y");
    titleName = sprintf("F%s",num2str(func_num));
    title(titleName);
end

