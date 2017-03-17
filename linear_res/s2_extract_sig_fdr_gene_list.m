clear
clc
fclose all;
mkdir res_report
letpaint=false;
regrlist=dir('regr*.txt');
cutoff=0.05;
ctf=num2str(cutoff);
fid3=fopen(['res_report/res_REGRESS_fdr',ctf(3:end),'_genelist.txt'],'w');
fid4=fopen(['res_report/res_REGRESS_fdr',ctf(3:end),'_summary.txt'],'w');

for tissueid=1:size(regrlist,1);
        disp(['now comes to ',regrlist(tissueid).name]);
        dirname=strrep(regrlist(tissueid).name,'regr_','');
        dirname=strrep(dirname,'.txt','');
        drimg=['res_report/',dirname,'/img'];
        mkdir(['res_report/',dirname,'/img']);
        fid1=fopen(['res_report/',dirname,'/res_REGRESS_fdr',ctf(3:end),'_genelist.txt'],'w');
        fid2=fopen(['res_report/',dirname,'/res_REGRESS_fdr',ctf(3:end),'_summary.txt'],'w');
        s_data_common;
        [geneidx,r,p,gene] = importfile([regrlist(tissueid).name]);
        % FDR correction
        q=mafdr(p,'BHFDR',true);
        
        % Get the residual values of the regression 
        residuals2=residuals';
        for k1=1:size(data,1)
            Y=data(k1,:)';
            W=Wxa(k1,:)';
            tmp=Y-covXa*W;
            residuals2(k1,:)=tmp';
        end
        Y=residuals2(geneidx,:)';
        [r2,p2]=corr(Y,PMI,'type','s');
        
        % FDR correction
        q2=mafdr(p2,'BHFDR',true);
        
        ix=q<0.05 & q2<cutoff;
        gid2=geneidx(ix);
        gen2=gene(ix);
        q1=q(ix);
        r1=r(ix);
        p2=p2(ix);
        r3=r2(ix);
        q3=q2(ix);
        k3=0;
        % Permutation test
        for k2=1:length(q3)
            Y=residuals2(gid2(k2),:)';
            for perm=1:10000
            perm_idx=randperm(length(PMI));
            PMIperm(:,perm)=PMI(perm_idx);
            end
            [~,pperm]=corr(PMIperm,Y,'type','s'); 
            clear PMIperm;
          % ------------------             
                if letpaint
                    Y=residuals2(gid2(k2),:)';
                    f=figure();
                    scatter(X,Y);
                    refline
                    box on
                    title(gen2{k2})
                    xlabel(sprintf('p=%.1e b=%.3f p=%.1e r=%.3f',q1(k2),r1(k2),q3(k2),r3(k2)));
                    saveas(f,['res_report/',tissuename,'/img/',tissuename,'_',ctf(3:end),'_',gen2{k2},'.jpg']);
                    close(f)
                end
	     % ------------------  
         % set the cutoff of permutations as 500 
            if sum(pperm<p2(k2))<500
                k3=k3+1;
            end
            fprintf(fid1,'%s\t%d\t%d\t%s\t%f\t%e\t%d\n',...
              tissuename,sum(ix),k2,gen2{k2},r3(k2),q3(k2),sum(pperm<p2(k2)));
            fprintf(fid3,'%s\t%d\t%d\t%s\t%f\t%e\t%d\n',tissuename,sum(ix),k2,gen2{k2},r3(k2),q3(k2),sum(pperm<p2(k2)));
        end     
     fclose(fid1);
     fprintf(fid2,'%s\t%d\t%d\t%d\t%d\t%d\t%d\t%d\n',tissuename,sum(q<0.05),sum(q2<cutoff),sum(r3>0),sum(r3<0),sum(ix),k3);
     fclose(fid2);
     fprintf(fid4,'%s\t%d\t%d\t%d\t%d\n',tissuename,sum(r3>0),sum(r3<0),sum(ix),k3);
     save(sprintf(['res_report/',dirname,'/%s_%s'],tissuename,ctf(3:end)),'gen2','r1','q1','r3','q3');
     disp(['Finished ',regrlist(tissueid).name]);
end

fclose(fid3);
fclose(fid4);
