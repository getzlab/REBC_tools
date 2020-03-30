function [counto,xco,yco,ho] = hist2(x,y,nx,ny,CLIM1,P)
% function [counto,xco,yco,ho] = hist2(x,y,nx,ny,CLIM1,P)
% HIST2 calculates a 2-dimensional histogram
%    [N,xb,yb,h] = HIST2(X,Y,NX,NY,CLIM,P) bins the X and Y data into bins NX and NY in
%    both dimensions. The color (z) scale is set by CLIM (optional) and a
%    colorbar is displayed unless P.BAR==false (optional).  counto is the count in each bin.  
%
%    N = HIST2(X,Y,M), where M is a scalar, bins the data into M equally
%    spaced bins in both dimensions
%
%    N = HIST2(X,Y,B), where B is a vector, bins the data with centers
%    specified by B 
%
%    The number of bins or centers can be specified individually for either
%    dimension using N = HIST2(X,Y,NX,NY) or N = HIST2(X,Y,BX,BY)
%
%    [N,BX,BY] = HIST2(...) also returns the positions of the bin centers
%    as two matrices in BX and BY
%
%    HIST2(...) without output arguments produces a colormapped image plot
%    of the 2d histogram
%
%    HIST2(...) also returns the positions of the bin centers
%    as two matrices in BX and BY
%
%    The last input argument is a struct of options. 
%       P.BAR=true (default) will make the colorbar 
%       P.MARGINHISTS=true (default)will make the margin histos.
%       P.COLORMAP= (default blue to red) 
%
%    The last output argument is a vector of handles to
%      ho(1) = 2d imagesc
%      ho(2) = colorbar
%      ho(3) = vertical histo on left
%      ho(4) = horizontal histo on bottom
%
% EXAMPLE
%   yin = randn(1,1000);
%   xin = randn(1,1000);
%   [n,x,y] = hist2d(xin,yin,11);
%   imagesc(x(1,:),y(:,1),n); hold on; plot(xin,yin,'y.'); colorbar

if ~exist('nx')
   nx = 10;
end

if ~exist('ny')
   ny = nx;
end

if ~exist('CLIM1')
   CLIM1 = NaN*[1 100];
end

if nargin<6
    P.BAR=true;
    P.MARGINHISTS=true;
end

if ~isfield(P,'BAR')
    P.BAR=true;
end    
if ~isfield(P,'MARGINHISTS')
    P.MARGINHISTS=true;
end
if ~isfield(P,'HIST2F')
    P.HIST2F=false;
end    
if ~isfield(P,'CUMULATIVE')
    P.CUMULATIVE=false;
end    
if ~isfield(P,'BOXPLOT')
    P.BOXPLOT=false;
end    

if ~isfield(P,'XLAB')
    P.XLAB='';
end    
if ~isfield(P,'YLAB')
    P.YLAB='';
end    
if ~isfield(P,'TITLE')
    P.TITLE='';
end    
if ~isfield(P,'HIST2F')
    P.HIST2F=false;
end    
if ~isfield(P,'BKG')
    P.BKG=1;
end   
if ~isfield(P,'COLORMAP')
    %v=(0:254)'/254;
    %cmap=[v 0*v 1-v]; cmap=[[1 1 1]; cmap];
    
    cmap=jet(1023);
    cmap(1,:)=P.BKG*[1 1 1];
    
    
    P.COLORMAP=cmap;
    
    if (P.HIST2F)
        P.COLORMAP=gray(254);
    end
end


    
if length(x) ~= length(y)
   error(sprintf('x and y must be same size ( %g ~= %g )',length(x),length(y)));
end
   
[dummy,xc] = hist(x,nx);
[dummy,yc] = hist(y,ny);

count = [];

for i = 1:length(yc)
   if i == 1
      lbound = -Inf;
   else
      lbound = (yc(i-1) + yc(i)) / 2;
   end
   if i == length(yc)
      ubound = inf;
   else
      ubound = (yc(i) + yc(i+1)) /2;
   end
   count(i,:) = hist(x((y >= lbound) & (y < ubound)),xc);
end

[xc, yc] = meshgrid(xc, yc);

hh=NaN;

if (nargout == 0)|| (nargout>3)
    %cmap=1-(1-jet).*((gray).^0.25);
    zc=count; zc(zc<1)=0.1;
    if any(isnan(CLIM1))|isempty(CLIM1)
        CLIM1=[0.9 10^(round(log10(max(zc(:))))+0)];
    end
    if (CLIM1(1)==1)
        CLIM1(1)=0.9;
    end
    if (CLIM1(2)==1)
        CLIM1(2)=2;
    end
    if (P.MARGINHISTS)
        subplot('position',[0.25 0.28 0.65 0.65])
    end
    h1=imagesc(xc(1,:),yc(:,1),log10(zc),log10(CLIM1));
    if (P.HIST2F)
        CLIM=[0 1]
        [nx1 ny1]=size(zc)
        for i=1:nx1
            zn(i,:)=zc(i,:)/sum(zc(i,:));
        end
        h1=imagesc(xc(1,:),yc(:,1),zn,[0 1]);        
    end
    
    set(gca,'ydir','normal')
    
    dx=diff(xc(1,1:2)); xlim1=[xc(1,1)-dx/2 xc(1,end)+dx/2];
    dy=diff(yc(1:2,1)); ylim1=[yc(1,1)-dy/2 yc(end,1)+dy/2];
        
     
    if P.CUMULATIVE
        [~,k]=sort(x(:));
        x=x(k);
        y=y(k);
        x1=NaN*y;   
        [~,b]=histc(x,[0 xc(1,:)]);
        u=unique(b);
        Nu=length(u);
        for i=1:Nu
            k=find(ismember(b,u(i)));
            [~,k1]=sort(y(k));
            y(k)=y(k(k1));
            x(k)=x(k(k1));
            nk=length(k);
            x1(k,1)=i-0.5+.95*(1:nk)/nk;    
        end        
        plot(x1,y,'o','markerfacecolor','b','markersize',5);
        set(gca,'xticklabel',[],'yticklabel',[])
        xL=repmat((1:Nu),2,1)+0.5; yL=repmat(ylim1',1,Nu);
        line(xL,yL,'linestyle','--','color',0.5*[1 1 1])        
        hh(1)=gca;
        hh(2)=NaN;
    elseif P.BOXPLOT
        [~,k]=sort(x(:));
        x=x(k);
        y=y(k);
        xb=[0 xc(1,:)]; xb(end)=xb(end)+1e-5;
        [~,b]=histc(x,xb);
        xu=unique(sort(b));
        Nu=length(xu);
        %subplot('position',[0.25 0.28 0.65 0.65])
        ax1=subplot('position',[0.305 0.317 0.6 0.61])
        boxplot(ax1,y,b,'notch','on')
        %set(gca,'XTickLabel',{' '})
        dy=diff(yc(1:2,1)); ylim1=[yc(1,1)-dy/2 yc(end,1)+dy/2];
        ylim(ylim1)
        p = kruskalwallis(y,x,'off')
        set(gca,'XTickLabel',{' '},'yticklabel',[])
        hh(1)=gca;
        hh(2)=NaN;        
    else
        hh(1)=get(h1,'parent');    
        %axis xy; 
        colormap(P.COLORMAP); %grid on;
        if (P.BAR)
            hh(2)=colorbar;
            if (~P.HIST2F)
                zt=get(hh(2),'ytick');
                lzt=round(10.^zt);
                zt=log10(unique(lzt));
                set(hh(2),'ytick',zt,'yticklabel',num2str(round(10.^zt')));
            end
        end
    end
    
    title(P.TITLE)
    
    if (P.MARGINHISTS)
        set(hh(1),'xticklabel',{' '},'yticklabel',[])
        subplot('position',[0.11 0.28 0.12 0.65])
        h1(3)=barh(yc(:,1),sum(count,2),'hist');
        
        hh(3)=get(h1(3),'parent');
        
        if ~(P.CUMULATIVE|P.BOXPLOT)
           set(gca,'ydir','reverse')
        end
        ylim(ylim1)
        set(hh(1),'ylim',ylim1)
        p1=get(hh(1),'position');
        ylabel(P.YLAB)
        if isfield(P,'YLABELS')
            set(gca,'ytick',1:length(yc(:,1)),'yticklabel',P.YLABELS)
        end
        
        set(hh(3),'ydir','normal')

        subplot('position',[0.25 0.117 p1(3) 0.12])
        h1(4)=bar(xc(1,:),sum(count,1),1,'EdgeColor',0.00005*[1 1 1]);%,'hist');
        hh(4)=get(h1(4),'parent');
        xlim(xlim1)
        set(hh(1),'xlim',xlim1)
        set(h1(3:4),'FaceColor',0.5*[1 1 1],'EdgeColor',0.25*[1 1 1])
        xlabel(P.XLAB)      
        if isfield(P,'XLABELS')
            set(gca,'xtick',1:length(yc(:,1)),'xticklabel',P.XLABELS)
        end
        
        linkaxes([hh(4) hh(1) ],'x');
        linkaxes([hh(3) hh(1) ],'y');

    end
end

if ~isnan(hh(1))
    axes(hh(1))
end

if (nargout>0)
    counto = count;
    xco = xc;
    yco = yc;
    ho=hh;
end

%%
function test
x=randn(5000,1);
y=randn(5000,1);
hist2D(x,y,-3:0.2:3,-3:0.2:3)
