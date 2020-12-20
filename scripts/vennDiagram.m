function vennDiagram(sv,lab,P) 
%function vennDiagram(sv,lab) 
%   plots a Venn based on input struct sv and lab

if (nargin<3)
    P.fs=14;
end
if ~isfield(P,'fs')
    P.fs=14;
end
fs=P.fs;

[nx ny]=size(sv.data);
R1=.97;
[X,Y] = cylinder(R1,1000);
if (ny==2)
    nlab(1) = cellstr(num2str(sum(sv.count(1+[1 3]))));
    nlab(2) = cellstr(num2str(sum(sv.count(1+[2 3]))));
end
if(ny==3)
    nlab(1) = cellstr(num2str(sum(sv.count(1+[1 3 5 7]))));
    nlab(2) = cellstr(num2str(sum(sv.count(1+[2 3 6 7]))));
    nlab(3) = cellstr(num2str(sum(sv.count(1+[4:7]))));
end
if(ny==5)
    nlab(1) = cellstr(num2str(sum(sv.count(1+[1:2:15]))));
    nlab(2) = cellstr(num2str(sum(sv.count(1+[2 3 6 7]))));
    nlab(3) = cellstr(num2str(sum(sv.count(1+[4:7]))));
end
nlab=strcat(lab,':',nlab);
text(1.6*R1,0,'.','fontsize',fs,'HorizontalAlignment','center')
text(-1.6*R1,0,'.','fontsize',fs,'HorizontalAlignment','center')
switch ny
    case 2
        line(X(1,:)-R1/2,Y(1,:),'linewidth',2)
        line(X(1,:)+R1/2,Y(1,:),'linewidth',2)
        %text(-R1/2,R1,nlab(1))
        %text(+R1/2,R1,nlab(2))
        text(-0.75*R1,1.1*R1,nlab(1),'fontsize',fs,'HorizontalAlignment','center')
        text(+0.75*R1,1.1*R1,nlab(2),'fontsize',fs,'HorizontalAlignment','center')
        set(gca,'visible','off')
        text(1.3*R1,-1.5*R1,num2str(sv.count(1)),'fontsize',fs,'HorizontalAlignment','center')
        text(-0.9*R1,0,num2str(sv.count(2)),'fontsize',fs,'HorizontalAlignment','center')
        text(+0.9*R1,0,num2str(sv.count(3)),'fontsize',fs,'HorizontalAlignment','center')
        text(0,0,num2str(sv.count(4)),'fontsize',fs,'HorizontalAlignment','center')
 
    case 3
        line(X(1,:)-R1/2,Y(1,:),'linewidth',2)
        line(X(1,:)+R1/2,Y(1,:),'linewidth',2)
        line(X(1,:),Y(1,:)-R1,'linewidth',2)
        text(-0.9*R1,1.1*R1,nlab(1),'fontsize',fs,'HorizontalAlignment','center')
        text(+0.9*R1,1.1*R1,nlab(2),'fontsize',fs,'HorizontalAlignment','center')
        text(0,-2.1*R1,nlab(3),'fontsize',fs,'HorizontalAlignment','center')
        set(gca,'visible','off')
        text(1.25*R1,-2.0*R1,num2str(sv.count(1)),'fontsize',fs,'HorizontalAlignment','center')
        text(-0.9*R1,0.4*R1,num2str(sv.count(2)),'fontsize',fs,'HorizontalAlignment','center')
        text(+0.9*R1,0.4*R1,num2str(sv.count(3)),'fontsize',fs,'HorizontalAlignment','center')
        text(0,0.3*R1,num2str(sv.count(4)),'fontsize',fs,'HorizontalAlignment','center')
        text(0,-1.5*R1,num2str(sv.count(5)),'fontsize',fs,'HorizontalAlignment','center')
        text(-0.5*R1,-0.7*R1,num2str(sv.count(6)),'fontsize',fs,'HorizontalAlignment','center')
        text(+0.5*R1,-0.7*R1,num2str(sv.count(7)),'fontsize',fs,'HorizontalAlignment','center')
        text(0,-0.35*R1,num2str(sv.count(8)),'fontsize',fs,'HorizontalAlignment','center')
                
end

%A = [300 300 300]; I = [100 100 100 100];
%venn(A,I,'FaceColor',{'w','w','w'},'EdgeColor','black')
 
end
    
function test
f='/Users/stewart/Projects/Cancer/UCS/Fusions/Algs.30Nov2014.txt'
F=load_tsv(f)
sv = vennCounts([F.aMDA F.aUNC F.aBI])
clf
vennDiagram(sv,{'MDA' 'UNC' 'BI'})
end