function [corrspec,corrtime,corrgood,label_list] = tfcomp1rotate_correctevent(TF_cal,TFs,transprop,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ);

TF_name = TF_cal{1};
TFidx = strmatch(TF_name,{TFs.label},'exact');
comp1 = TF_name(1);
comp2 = TF_name(2);
label_list{1} = [comp1,comp2];

if strcmp(comp1,'Z')==1
    spec_1= spec_mat(:,1);
    isgood1 = goodZ;
else
    disp('Error, invalid component for rotation')
    return
end

if strcmp(comp2,'H')==1
    spec_H1 = spec_mat(:,2);
    spec_H2 = spec_mat(:,3);
    isgood2 = goodH;
else
    disp('Error, invalid component for rotation')
    return
end

hang = transprop.params.hang;
cang = cos(hang*pi/180);
sang = sin(hang*pi/180);

spec_2 = sang.*spec_H2+cang.*spec_H1;

corrspec = spec_1-((TFs(TFidx).transfunc_tap)'.*spec_2);
amp1_2 = real(ifft(2*corrspec,NFFT)./dt);
corrtime = amp1_2(npad0+1:npad0+npts);

corrgood = isgood1*isgood2;

return