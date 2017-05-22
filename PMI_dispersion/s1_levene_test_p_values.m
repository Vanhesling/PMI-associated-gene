clc
clear
fclose all;
mkdir res

for tissueid=1
%% Define the S-PMI and L-PMI groups     
    s_data_common;
    isPMI=PMI>=median(PMI);
    n=size(PMI,1);
    [~,i]=sort(PMI);
    PMIgroup=2*ones(size(PMI));
    PMIgroup(i(1:round(n/3)))=1;
    PMIgroup(i(end-round(n/3)+1:end))=3;
%% Levene's test

fid=fopen(sprintf('res/levene_res_%s.txt',tissuename),'w');

for k=1:length(g_id)  
    y1=data(g_id(k),:)';
    W=Wxa(g_id(k),:)';
    Y=y1-covXa*W;
    p1=vartestn(Y,PMIgroup,'TestType','LeveneAbsolute','display','off');
    p2=vartestn(Y,PMIgroup,'TestType','LeveneQuadratic','display','off');
    p3=vartestn(Y,PMIgroup,'TestType','BrownForsythe','display','off');
    fprintf(fid,'%d\t%e\t%e\t%e\t%d\t%s\n',g_id(k),p1,p2,p3,...
        double(var(Y(PMIgroup==1))>var(Y(PMIgroup==3))),g{g_id(k)});
end

fclose(fid);
end