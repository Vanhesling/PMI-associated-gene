
clear
clc
fclose all;
load('E:\My_project\3.results\peer_exclude_noneuropean\Genotype PMI interaction\list450idv');

tissueid=18;
s_data_common;
[newid,idx1,idx2]=intersect(uid,id);
modelstr=['y~1+',sprintf('x%d+',1:size(covXb,2)+2),'x1:x2'];
residuals2=data'-(covXa*Wxa');
residuals2=residuals2(idx1,:);
pmi=PMI(idx1);
covXb=covXb(idx1,:);
%% Genotype-PMI-Interaction

% tissueid=18; chrid=14; targetgene='NAP1L2'; snppos=80165382; alleles1='GG',alleles2='GA/AA',alleles3='GA',alleles4='AA',rsidx='rs11159412';
% tissueid=18;chrid=1;targetgene='SPIN3',snppos=151782130;alleles1='GG';alleles2='GT/TT',alleles3='GT',alleles4='TT',rsidx='rs1521177';

for chrid=1:22
    chrid
    load(['E:\My_project\3.results\peer_exclude_noneuropean\Genotype PMI interaction\mat\GTEx_20150112_450Indiv_chr',num2str(chrid)]);
    load(['E:\My_project\3.results\peer_exclude_noneuropean\Genotype PMI interaction\mat\mark_GTEx_20150112_450Indiv_chr',num2str(chrid)]);    
    geno=geno(idx2,:);
    mafv=snp_maf(geno);
    [geno,mark]=snp_pickmarker(geno,mark,mafv>0.15);
    g012=snp_012geno(geno);
    [~,yid]=intersect(genelist,targetgene);
    y=residuals2(:,yid); 
        for gk=1:size(g012,2)
            gk=find(mark.pos==snppos);
            snp=g012(:,gk);
            mdl = fitlm([pmi,snp,covXb],y,modelstr);
            interactive_p=mdl.Coefficients.pValue(end);
            if interactive_p<1e-7
               fprintf(fid,'%s\t%s\t%d\t%d\t%s\t%e\n',...
               tissuename,mark.rsid{gk},chrid,mark.pos(gk),genelist{yid},interactive_p);
            end
        end
end
     
%% Plot the PMI-Genetype Interaction
 
        color1='[0 0 .7]';
        color2='[0.7 0 0]';
        close all
        h=figure;
        subplot(1,4,1)
        sc1=scatter(pmi(snp==0),y(snp==0),'filled');
        hold on
        sc2=scatter(pmi(snp>=1),y(snp>=1),'filled');

        xlabel('PMI (min)','FontSize',12,'FontWeight','bold')
        ylabel({targetgene;'residual expression'},'FontSize',14,'FontWeight','bold')
        sc1.MarkerEdgeColor=color1;
        sc1.MarkerFaceColor=color1;
        sc1.MarkerFaceAlpha=0.5;
        sc2.MarkerEdgeColor=color2;
        sc2.MarkerFaceColor=color2;
        sc2.MarkerFaceAlpha=0.5;
        %legend({alleles1,alleles2},'location','northoutside','Orientation','horizontal')

         rf=refline;
         rf(1).Color=color2;
         rf(1).LineStyle=':';  
         rf(2).Color=color1;
         rf(2).LineStyle=':';

        box on
        xlim([0 1600])
        plot(pmi(snp==2),y(snp==2),'kx')
        yy=[-0.5 0.5];
        ylim(yy);
        R=corr(pmi,y,'type','p');
        text(600,-0.4, ['Pearson-R = ' num2str(R)])
        title(rsidx,'FontSize',16,'FontWeight','bold')


        subplot(1,4,2)
        sc1=scatter(pmi(snp==0),y(snp==0),'filled'); box on; 
        R=corr(pmi(snp==0),y(snp==0),'type','p');
        text(600,-0.4, ['Pearson-R = ' num2str(R)])
        legend(sprintf('%s - %s',alleles1,rsidx),'location','northoutside','Orientation','horizontal')
        xlim([0 1600])
        ylim(yy);
        rf=refline;
        rf.Color=color1;
        rf.LineStyle=':';
        xlabel('PMI (min)','FontSize',12,'FontWeight','bold')
        sc1.MarkerEdgeColor=color1;
        sc1.MarkerFaceColor=color1;
        sc1.MarkerFaceAlpha=0.5;

        subplot(1,4,3)
        sc2=scatter(pmi(snp==1),y(snp==1),'filled'); box on; refline
        sc2.MarkerEdgeColor=color2;
        sc2.MarkerFaceColor=color2;
        sc2.MarkerFaceAlpha=0.5;
        R=corr(pmi(snp==1),y(snp==1),'type','p');
        text(600,-0.4, ['Pearson-R = ' num2str(R)])
        legend(sprintf('%s - %s',alleles3,rsidx),'location','northoutside','Orientation','horizontal')
        xlim([0 1600])
        ylim(yy);
        rf=refline;
        rf.Color=color2;
        rf.LineStyle=':';
        xlabel('PMI (min)','FontSize',12,'FontWeight','bold')

        
        subplot(1,4,4)
        sc2=scatter(pmi(snp==2),y(snp==2),'filled'); box on; refline
        sc2.MarkerEdgeColor=color2;
        sc2.MarkerFaceColor=color2;
        sc2.MarkerFaceAlpha=0.5;
        R=corr(pmi(snp==2),y(snp==2),'type','p');
        text(600,-0.4, ['Pearson-R = ' num2str(R)])
        legend(sprintf('%s - %s',alleles4,rsidx),'location','northoutside','Orientation','horizontal')
        xlim([0 1600])
        ylim(yy);
        rf=refline;
        rf.Color=color2;
        rf.LineStyle=':';
        xlabel('PMI (min)','FontSize',12,'FontWeight','bold')
        hold on
        plot(pmi(snp==2),y(snp==2),'kx')
        %%
        %set(h, 'Position', [0, 0, 350, 1100]);
        set(h, 'PaperPosition', [0.635, 6.35, 40.32, 8.24]);
        %set(gcf, 'PaperPosition', [2 1 10 2]);
        print(h,targetgene, '-dpng', '-r1200');
 
     
     

