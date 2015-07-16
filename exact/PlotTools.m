function h = PlotTools(  )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

h.OverlayPColor = @OverlayPColor;
h.SurfaceDiffPlot = @SurfaceDiffPlot;





end


function SurfaceDiffPlot(grid,l1,l2)

    e=grid.E(:,l2)-grid.E(:,l1);
    dim=sqrt(length(e));
    re=reshape(e,dim,dim);
    surf(re);
    

end

function [Zi c2 c3] = OverlayPColor(grid,con,sh,n)

dim=[length(unique(grid.x)),length(unique(grid.y))];

[x y th a1 a2] = ploch(grid,dim,n);

%x=gridshift(x,sh);
%y=gridshift(y,sh);
th=gridshift(th,sh);
a1=gridshift(a1,sh);
a2=gridshift(a2,sh);
 %= meshgrid(0:.01:2*pi,0:.01:2*pi);
 xmax = max(max(x));
 xmin = min(min(x));
 ymax = max(max(y));
 ymin = min(min(y));
[ Xi Yi] = meshgrid(xmin:.01:xmax,ymin:.01:ymax);


A1 = interp2(x,y,a1,Xi,Yi,'linear');
A2 = interp2(x,y,a2,Xi,Yi,'linear');

Zir = interp2(x,y,real(th),Xi,Yi);
Zii = interp2(x,y,imag(th),Xi,Yi);
Zi=angle(Zir+Zii*1i);

%Xi=x;
%Yi=y;
%Zi=th;
%Ai=a2;

hs=pcolor(Xi,Yi,Zi);
shading interp;
hsa=gca;
crange=caxis;
hold on
colormap hsv;
xlabel('\theta_x');
ylabel('\theta_y');
set(gca,'XTick',[0+.1,pi,2*pi-.001])
set(gca,'XTickLabel',{'-pi','0','pi'})
set(gca,'YTick',[0+.1,pi,2*pi-.001])
set(gca,'YTickLabel',{'-pi','0','pi'})


[c hc] = contour(x,y,log10(a1),[con,con],'w');
    hlines=findobj(hc,'Type','line');
    set(hc,'LineWidth',2);
    
%[c2 hc2] = contour(Xi,Yi,log10(A1),[con(1),con(1)],'k');
%set(hc2,'LineWidth',2);
%[c3 hc3] = contour(Xi,Yi,log10(A2),[con(2),con(2)],'r');
%set(hc3,'LineWidth',2); 
%caxis([-pi,pi]);
hold off


end

function [ x y th A1 A2 ] = ploch( ch ,dim,n)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

%dim=size(ch.x);

th=reshape(ch.th(:,n),dim(1),dim(2));
x=reshape(ch.x,dim(1),dim(2));
y=reshape(ch.y,dim(1),dim(2));
A1=abs(reshape(ch.A1(:,n),dim(1),dim(2)));
A2=abs(reshape(ch.A2(:,n),dim(1),dim(2)));

%surf(x,y,log10(A2),th,'LineStyle','none');
%colormap hsv
%xlabel('\theta_x')
%ylabel('\theta_y')
end


function X = gridshift(X,dx)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

d=size(X);

if dx(1) < 0
    dx(1)=abs(dx(1));
    X=[X(:,(dx(1)+1):d)';X(:,1:dx(1))']';

elseif dx(1) > 0
    X=[X(:,1:dx(1))';X(:,(dx(1)+1):d)']';    
end 
    

if dx(2) < 0
    dx(2)=abs(dx(2));
    X=[X((dx(2)+1):d,:);X(1:dx(2),:)];

elseif dx(2) > 0
    X=[X(1:dx(2),:);X((dx(2)+1):d,:)];    
end 






end

function [Xi Yi] = XYINTERP(X,Y,res)
%UNTITLED2 take an array and fill in points to increase resolution.
%   takes arguments (fac,per, th, x, y). th is the function, x and y are
%   coordiants, and are optional. fac in the number of intersital points to
%   insert. pbc is a logical, true if the function th has beriodic boundary
%   ocnditions

xii=(min(min(X)):res:max(max(X)));
nxd=length(xii);
yii=(min(min(Y)):res:max(max(Y)));
nyd=length(yii);

Xi=ones(nyd,1)*xii;
Yi=yii'*ones(1,nxd);

end



