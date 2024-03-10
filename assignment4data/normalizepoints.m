function [xn,N]=normalizepoints(x);
%input homogenous points x(:,i)=[xi;yi;1]
%output normalized homogenous points xn
% and normalization matrix N

%calculate mean and std of points
xmean=mean(x(1:2,:),2); %(x0,y0)
xstd=std(x(1:2,:),0,2); %(sx,sy)

%Normalization matrix
N=[1/xstd(1), 0, -xmean(1)/xstd(1);0, 1/xstd(2), -xmean(2)/xstd(2);0, 0, 1];

if(size(x,1)==2)
    y=[x;ones(1,length(x))];
else
    y=x;
end

xn=N*y;
end