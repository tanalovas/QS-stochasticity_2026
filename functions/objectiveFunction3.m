function errorFunction = objectiveFunction3(X,data,flag,flagPlot)

toll = 0.1e-1;
X    = abs(X);

% Data to be fitted 
%   1      2       3         4       5      6        7       8       9      10
% time C8_mean HC10_mean HC8_mean C8_min C8_max HC10_min HC10_max HC8_min HC8_max
tt      = data(:,1);
c8      = data(:,2);
hc10    = data(:,3);
hc8     = data(:,4);
c8neg   = data(:,5);
c8pos   = data(:,6);
hc10neg = data(:,7);
hc10pos = data(:,8);
hc8neg  = data(:,9);
hc8pos  = data(:,10);


[tsol,ysol] = solveODEfit(X,data);

index = zeros(size(tt));

for i = 1:length(index)
    for j = 1:length(tsol)
        if(abs(tsol(j)-tt(i)) < toll)
            index(i) = j;
            break
        end
    end
end


ysolMini = ysol(index,:);

errorC8   = (1/length(tt))*sum((ysolMini(:,end)-c8).^2);
errorHC10 = (1/length(tt))*sum((ysolMini(:,end-1)-hc10).^2);
errorHC8  = (1/length(tt))*sum((ysolMini(:,end-2)-hc8).^2);

switch flag
    case 0
        errorFunction = (10*errorC8+errorHC10+errorHC8)/3;
    case 1
        errorFunction = errorC8;
    case 2
        errorFunction = errorHC10;
    case 3
        errorFunction = errorHC8;
end

if(flagPlot)
    figure(3);
    subplot(3,1,1)
    plot(tsol,ysol(:,end),'k',tsol(index),ysolMini(:,end),'k*',tt,c8,'*r',tt,c8neg,'xg',tt,c8pos,'xg');
    subplot(3,1,2)
    plot(tsol,ysol(:,end-1),'k',tsol(index),ysolMini(:,end-1),'k*',tt,hc10,'*r',tt,hc10neg,'xg',tt,hc10pos,'xg');
    subplot(3,1,3)
    plot(tsol,ysol(:,end-2),'k',tsol(index),ysolMini(:,end-2),'k*',tt,hc8,'*r',tt,hc8neg,'xg',tt,hc8pos,'xg');
end

