function [h,sv] =  upset_plot(X,algs,P)

if (nargin<3)
    P=[]
    P.log=0
    P.all=0
    P.angle=45
end
if ~isfield(P,'log')
    P.log=0
end
if ~isfield(P,'all')
    P.all=0
end
if ~isfield(P,'angle')
    P.angle=45
end

NA=length(algs)
q=[]
for i=1:length(algs)
    q(:,i)=X.(algs{i})>0;
end
sv = vennCounts(q)


% clf
[n k]=sort(sv.count,'descend')
b=(1:length(sv.count))-1; b=b(k)
NB=length(b)
h(1) = axes( 'Position', [0.2 0.3 0.6 0.6] );
hb=bar(h(1),n);
set(h(1),'xticklabels',[])
grid on
xlim(0.5+[0 NB])
la=n(1)/50;
lb=1.0;
if (P.log)
    set(h(1),'yscale','log')
    set(hb,'BaseValue',0.9)
    ylim([0.9 max(n)*1.5])
    la=0;
    lb=1.2;
else
    ylim([0 max(n)*1.05])    
end

ylabel('events')

NB1=NB;
if (P.all<1)
    NB1=sum(n>0)
    xlim(0.5+[0 NB1])
    set(h(1),'xticklabels',[])    
end

for i=1:NB1
    text(i,la+n(i)*lb,sprintf('%d',n(i)),'fontsize',9,'horizontalalignment','left','rotation',P.angle,'color','k')
end



h(2) = axes( 'Position', [0.2 0.15 0.6 0.14],'XColor', [1,1,1], 'YColor', 1*[1,1,1] );
axis([0.5+[0 length(sv.count)] 0.5+[0 NA]])
xlim(0.5+[0 NB1])
grid on
%set(h(2),'xticklabels',[])% ,'yticklabels',[])%,'visible','off')
set(h(2),'ytick',1:NA,'yticklabels',algs)
text(-0.05+zeros(NA,1),[1:NA],algs,'horizontalalignment','right','fontsize',14)
hold on
for i=1:NB1
    b1=dec2bin(b(i))
    k=strfind(fliplr(b1),'1')
    text(i+zeros(size(k)),k,'+','fontsize',14,'horizontalalignment','center')
    plot(i+zeros(size(k)),k,'o','markersize',10,'markerfacecolor','k','markeredgecolor','k')
    line(i+[0 0],[1 NA],'linestyle','-','color','k')
end
hold off
if (P.all<1)
    xlim(0.5+[0 NB1])
end

h(3) = axes( 'Position', [0.81 0.15 0.14 0.14] );
na=sum(sv.data)
hh=barh(h(3),na);
set(h(3),'yticklabels',[])
set(gca, 'XAxisLocation', 'top')
xlabel('events')
grid on
ylim([0.5+[0 NA]])
if (P.log)
    set(h(3),'xscale','log')
    set(hh,'BaseValue',0.9)
    xlim([0.9 max(na)*1.2])
    lx=ceil(log10(max(na)))
    set(h(3),'xtick',10.^[0:lx])
end
        
