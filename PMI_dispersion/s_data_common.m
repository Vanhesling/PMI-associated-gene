
filelist=dir('../expr*.mat');
tissuename={filelist.name};
tissuename=strrep(tissuename,'expr_','');
tissuename=strrep(tissuename,'.mat','');
tissuename=tissuename{tissueid};
tissuetextshort={'Adipose','Aorta','Tibial','Cerebellum','Cerebellum','Esophagus','Heart','Lung','Muscle','Nerve','Pituitary','Skin Suprapubic','Skin Lowerleg','Thyroid','Whole Blood'};

%% loadfile
    load(['..\expr_',tissuename]);
    load(['..\peer_ ',tissuename]);
    load(['..\smpl_',tissuename]);
    load('..\genelist.mat');
%% Define parameters
    g=genelist;
    PMI=double(uischtm);  
% filter the lowerly expressed genes (20%)
    isvalid_geneidx=(sum(data==0,2)/size(data,2))<0.5;
    mx=median(data,2);
    ishighly_exr=mx>prctile(mx,20);
    g_idx=ishighly_exr&isvalid_geneidx;
    g_id=[];
    for kk=1:length(g) 
        if g_idx(kk)
           g_id=[g_id kk];
        end
<<<<<<< HEAD
    end

% exclude the factors showing a Pearson¡¯s correlation or Spearman¡¯s rank correlation test P-value smaller than 0.05 
=======
    end 
% exclude the factors showing a Pearson's correlation or Spearman's rank correlation test P-value smaller than 0.05 
>>>>>>> origin/master
    covX=factors(:,6:end);
    Wx=weigthx(:,6:end);
    [~,p1]=corr(covX,PMI,'type','s');
    [~,p2]=corr(covX,PMI,'type','p');
    ix=p1<0.05|p2<0.05;
    covX(:,ix)=[];
    Wx(:,ix)=[];

% store the covariates and weights
    covXa=[factors(:,2:5) covX];
    Wxa=[weigthx(:,2:5) Wx];
    covXb=[factors(:,2:4) covX];
    Wxb=[weigthx(:,2:4) Wx];
    XX=[PMI covXb];
  
