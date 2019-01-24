function [corrspec1_2,corrspec1_32,corrspec1_432,corrtime1_2,corrtime1_32,corrtime1_432,corrgood1_2,corrgood1_32,corrgood1_432,label_list] =...
    tfcomp3_correctevent(TF_cal,TFs,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ)

TF_name = TF_cal{1};

comp1 = TF_name(1);
comp2 = TF_name(5);
comp3 = TF_name(4);
comp4 = TF_name(2);

TFidx12 = strmatch([comp1,comp2],{TFs.label},'exact');
TFidx32 = strmatch([comp3,comp2],{TFs.label},'exact');
TFidx42 = strmatch([comp4,comp2],{TFs.label},'exact');
TFidx43_2 = strmatch([comp4,comp3,'-',comp2],{TFs.label},'exact');
TFidx13_2 = strmatch([comp1,comp3,'-',comp2],{TFs.label},'exact');
TFidx14_32 = strmatch([comp1,comp4,'-',comp3,comp2],{TFs.label},'exact');

label_list{1} = [comp1,comp2];
label_list{2} = [comp1,comp3,'-',comp2];
label_list{3} = [comp1,comp4,'-',comp3,comp2];

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

if strcmp(comp3,'Z')==1
    spec_3 = spec_mat(:,1);
    isgood3 = goodZ;
elseif strcmp(comp3,'1')==1
    spec_3 = spec_mat(:,2);
    isgood3 = goodH;
elseif strcmp(comp3,'2')==1
    spec_3 = spec_mat(:,3);
    isgood3 = goodH;
elseif strcmp(comp3,'P')==1
    spec_3 = spec_mat(:,4);
    isgood3 = goodP;
end

if strcmp(comp4,'Z')==1
    spec_4 = spec_mat(:,1);
    isgood4 = goodZ;
elseif strcmp(comp4,'1')==1
    spec_4 = spec_mat(:,2);
    isgood4 = goodH;
elseif strcmp(comp4,'2')==1
    spec_4 = spec_mat(:,3);
    isgood4 = goodH;
elseif strcmp(comp4,'P')==1
    spec_4 = spec_mat(:,4);
    isgood4 = goodP;
end

corrspec1_2 = spec_1-((TFs(TFidx12).transfunc_tap)'.*spec_2);
% corrspec3_2 = spec_3-((TFs(TFidx32).transfunc_tap)'.*spec_1);
corrspec3_2 = spec_3-((TFs(TFidx32).transfunc_tap)'.*spec_2);
corrspec1_32 = corrspec1_2-((TFs(TFidx13_2).transfunc_tap)'.*corrspec3_2);
corrspec4_2 = spec_4-((TFs(TFidx42).transfunc_tap)'.*spec_2);
corrspec4_32 = corrspec4_2-((TFs(TFidx43_2).transfunc_tap)'.*corrspec3_2);
corrspec1_432 = corrspec1_32-((TFs(TFidx14_32).transfunc_tap)'.*corrspec4_32);

amp1_2 = real(ifft(2*corrspec1_2,NFFT)./dt);
amp1_32 = real(ifft(2*corrspec1_32,NFFT)./dt);
amp1_432 = real(ifft(2*corrspec1_432,NFFT)./dt);

corrtime1_2 = amp1_2(npad0+1:npad0+npts); 
corrtime1_32 = amp1_32(npad0+1:npad0+npts);
corrtime1_432 = amp1_432(npad0+1:npad0+npts); 

corrgood1_2 = isgood1*isgood2;
corrgood1_32 = isgood1*isgood2*isgood3;
corrgood1_432 = isgood1*isgood2*isgood3*isgood4;


return