clear
clc
fclose all;
filelist=dir('../expr*.mat');
for tissueid=1:length(filelist)
    s_data_common;
%% Multiple linear regression
    fid=fopen(['regr_',tissuename,'.txt'],'w');
    disp(['now comes to ',tissuename]);
    for k=1:length(g_id)
        Y=data(g_id(k),:)';
        pval=fitlm(XX,Y);
        pvalue=pval.Coefficients.pValue(2);
        estimate=pval.Coefficients.Estimate(2);
        fprintf(fid,'%d\t%f\t%e\t%s\n',g_id(k),estimate,pvalue,g{g_id(k)});
    end
    disp(['finish ',tissuename]);
    fclose(fid);
end





    
    
    
