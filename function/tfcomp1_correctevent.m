function [corrspec,corrtime,corrgood,label_list] = tfcomp1_correctevent(TF_cal,TFs,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ)

TF_name = TF_cal{1};
TFidx = strmatch(TF_name,{TFs.label},'exact');
comp1 = TF_name(1);
comp2 = TF_name(2);
label_list{1} = [comp1,comp2];

if strcmp(comp1,'Z')==1
    spec_1 = spec_mat(:,1);
    isgood1 = goodZ;
elseif strcmp(comp1,'1')==1
    spec_1 = spec_mat(:,2);
    isgood1 = goodH;
elseif strcmp(comp1,'2')==1
    spec_1 = spec_mat(:,3);
    isgood1 = goodH;
elseif strcmp(comp1,'P')==1
    spec_1 = spec_mat(:,4);
    isgood1 = goodP;
end

if strcmp(comp2,'Z')==1
    spec_2 = spec_mat(:,1);
    isgood2 = goodZ;
elseif strcmp(comp2,'1')==1
    spec_2 = spec_mat(:,2);
    isgood2 = goodH;
elseif strcmp(comp2,'2')==1
    spec_2 = spec_mat(:,3);
    isgood2 = goodH;
elseif strcmp(comp2,'P')==1
    spec_2 = spec_mat(:,4);
    isgood2 = goodP;
end

corrspec = spec_1-((TFs(TFidx).transfunc_tap)'.*spec_2); 
amp1_2 = real(ifft(2*corrspec,NFFT)./dt);
corrtime = amp1_2(npad0+1:npad0+npts);

corrgood = isgood1*isgood2;

return