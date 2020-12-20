function [h,xbin,Stats]=myhist(x,b,arg)
%MYHIST  Histogram plot.
%   [N,XBIN,STATS] = MYHIST(Y,B,ARG) adds info to histogram plot
%   input:
%   X,B: same as HIST
%   ARG = 'log' : log10 y scale
%         'gauss' : draws line with Gaussian fit
%         'nobox' : no stat box
%         'leftbox' : left stat box
%         'xlog'  : log x scale (remove *'s) 
%   output:
%   H,XBIN: same as hist
%   STATS = structure with mean std min max median ovr undr
%
%   Class support for inputs Y, X: 
%      float: double, single
%
%   See also HIST.

if (nargin<3)
    arg=false(1);
end

xlog=0;
if strfind(arg,'log')
    if strfind(arg,'xlog')
       xlog=1;
       arg=regexprep(arg,'xlog','')
       if length(b)<2
           b=logspace(round(min(log10(x)))-1,round(max(log10(x)))+1,b);
       end
    end
end
leftbox=0;
if strfind(arg,'left')
    leftbox=1;
end

if iscell(x)
    q=tab(x);
    h=bar(q.n);
    xbin=1:length(q.n);
    Stats.tot=sum(q.n);
    set(gca,'xtick',1:length(q.n),'xticklabel',q.x);
    return;
end

x=double(x);
[h,xbin]=hist(x,b);

if strfind(lower(arg),'cum')
  h=cumsum(h)
  if strfind(lower(arg),'norm')
    h=h/h(end);
  end
end

if (~isempty(arg))||(nargout==0)
    if(ishold)
        h01=findobj(gca,'Type','patch');
    else
        h01=[];
    end
    bar(xbin,h,minmax(xbin),'hist');
    h02=findobj(gca,'Type','patch');
    h0=setdiff(h02,h01);
    for h1=h0(:)'
        v0=get(h1,'Vertices'); v1=v0;
        if strfind(arg,'log')
            set(gca,'yscale','log');
            v1(v0(:,2)==0,2)=.8; set(h1,'Vertices',v1);
        end
        if (~isempty(h0))&&(ishold)
            %cc=length(h02)-1;
            co=get(0,'DefaultAxesColorOrder');
            c1=co(length(h01)+1,:);
            set(h1,'edgecolor',c1,'facecolor',c1);
        end
    end
   if (~ishold)
        if strfind(arg,'log')
          set(gca,'ylim',[max([0.8 min(h(:))])  10^(1.1*log10(max(h(:))))]);
        else
            set(gca,'ylim',[max([0 min(h(:))])  1.1*max(h(:))]);
        end
        if strfind(arg,'gauss')
            xlim=minmax(xbin);
            xx=x( (x>=xlim(1))&(x<=xlim(2)) );
            mu=mean(xx);
            sigma=std(double(xx));
            dx=diff(xbin(1:2));
            g= normpdf(xbin,mu,sigma);
            g=length(xx)*g/sum(g);
            hold on; plot(xbin,g,'k--'); hold off;
        end
    end
end

if (xlog)
    set(gca,'xscale','log')
    hL = findobj(gca,'Type','line');
    set(hL,'Visible','off') 
end


s=sprintf('N=%d',length(x));
xslim=minmax(double(x(:)'));
if strfind(arg,'statlim')
    xslim=minmax(xbin);
end
Stats.tot=length(x( (x>=xslim(1))&(x<=xslim(2)) ) );
Stats.mean=mean(x( (x>=xslim(1))&(x<=xslim(2)) ));
Stats.median=median(x( (x>=xslim(1))&(x<=xslim(2)) ));
Stats.std=std(double(x( (x>=xslim(1))&(x<=xslim(2)) )));
s=strvcat(s,sprintf('\\mu=%s',num2str(Stats.mean,' %.4g ')));
s=strvcat(s,sprintf('\\sigma=%s',num2str(Stats.std,' %.4g ')));
if length(b)>1
    %Stats.mean=mean(x( (x>=b(1))&(x<=b(end))));
    %Stats.std=std(double(x( (x>=b(1))&(x<=b(end)))));
    Stats.over=sum( x>b(end));
    Stats.under=sum( x<b(1));
    s=strvcat(s,sprintf('ovr=%s',num2str(Stats.over,' %d ')));
    s=strvcat(s,sprintf('und=%s ',num2str(Stats.under,' %d ')));
    db=b(2)-b(1);
    set(gca,'xlim',[b(1)-db  b(end)+db]);
else
    Stats.over=0;
    Stats.under=0;
end
%if ((~isempty(arg))||(nargout==0))&&(~ishold)&&((~isnumeric(arg))&&(~isempty(strfind(arg,'nobox'))))
if (((~isempty(arg))||(nargout==0))&&(~ishold))&&(isempty(strfind(arg,'nobox')))

    xb=0.97; 
    ha='right';
    if leftbox>0
        xb=0.03;
        ha='left';
    end
    h=text('units','normalized','position',[xb 0.97],'string',cellstr(s),...
    'EdgeColor','b','HorizontalAlignment',ha,'VerticalAlignment','top','interpreter','tex','backgroundcolor','w');
    grid on;

end
if nargout == 0
  clear h X* St*
end

