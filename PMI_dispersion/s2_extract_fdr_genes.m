clear 
clc
fclose all;
addpath('res');
a=dir('res/levene_res_*.txt');
mkdir res_report
mkdir img
fid=fopen('res_report/GENE_fdr05x_res_summary.txt','w');

for k=1:size(a,1) 
    [geneidx,p,p1,p2,r,gene] = importfile(a(k).name);
%   FDR correction
    q=mafdr(p,'BHFDR',true);
    ix=q<0.05&p2<0.05;
    if ~isempty(ix)
        id2=geneidx(ix);
        r2=r(ix);
        p2=p(ix);
        gene2=gene(ix);
        fprintf(fid,'%d\t%d\t%d\t%s\n',length(id2),sum(r2==1),sum(r2==0),a(k).name);
    else
        fprintf(fid,'%d\t%d\t%d\t%s\n',0,0,0,a(k).name);
    end
end
fclose(fid);

% FDR correction and Plot the DV genes
a=dir('res/levene_res_*.txt');
fid=fopen('res_report/GENE_fdr05x_res_list.txt','w');
for k=1:size(a,1) 
    tissueid=k;  
    s_data_common;
    [geneidx,p,p1,p2,r,gene] = importfile(a(k).name);
    [p,i]=sort(p);
    geneidx=geneidx(i);
    p2=p2(i);r=r(i);gene=gene(i);q=mafdr(p,'BHFDR',true);
    
    ix=q<0.05&p2<0.05;
    p=p(ix);p2=p2(ix);r=r(ix);
    gene=gene(ix);geneidx=geneidx(ix);
    q=q(ix);    
    
    isLPMI=PMI>=median(PMI);
    n=size(PMI,1);
    [~,ii]=sort(PMI);
    PMIgroup=2*ones(size(PMI));
    PMIgroup(ii(1:round(n/3)))=1;
    PMIgroup(ii(end-round(n/3)+1:end))=3;
    if ~isempty(ix)
        for j=1:length(q)
            y1=data(geneidx(j),:)';
            W=Wx(geneidx(j),:)';
            y=y1-covX*W;       
            
%             hh=figure('Visible','off');
%             subplot(1,2,1)
%             isLPMI=PMI>=median(PMI);
%             boxplot(y,isLPMI,'symbol','r+','colors',[0 0.3 0.5]);
%             hold on
%             isLPMI2=isLPMI+(rand(size(isLPMI))-0.5)*0.1+1;    
%             plot(isLPMI2,y,'o')
%             ss1=sprintf('%d-%d',min(PMI(~isLPMI)),max(PMI(~isLPMI)));
%             ss2=sprintf('%d-%d',min(PMI(isLPMI)),max(PMI(isLPMI)));
%             set(gca,'XTickLabel',{ss1,ss2})
%             xlabel('PMI group (min)','FontSize',12,'FontWeight','bold')
%             ylabel({gene{j};'residual expression'},'FontSize',12,'FontWeight','bold')
%             title(tissuetextshort{tissueid},'FontSize',16,'FontWeight','bold')

%             subplot(1,2,2)
%             plot(PMI,y,'o')
%             hold on;
%             plot([median(PMI) median(PMI)],get(gca, 'YLim'), 'r:','LineWidth', 1) 
%             xlabel('PMI (min)','FontSize',12,'FontWeight','bold')
%             ylabel({gene{j};'residual expression'},'FontSize',12,'FontWeight','bold')
%             title(tissuetextshort{tissueid},'FontSize',16,'FontWeight','bold')
%           
%             fname=sprintf('img/%s_%d_%d_%s',tissuename,j,r(j),gene{j});        
%             set(hh, 'PaperPosition', [0.635, 6.35, 25.32, 12.24]);
%             print(hh,fname,'-dpng','-r1200');
%             close(hh);          
            fprintf(fid,'%s\t%s\t%e\t%e\t%d\n',a(k).name,gene{j},p(j),q(j),r(j));
        end
    end    
end
fclose(fid);




