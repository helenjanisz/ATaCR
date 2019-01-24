function [corrspec1_2,corrspec1_32,corrtime1_2,corrtime1_32,corrgood1_2,corrgood1_32,label_list] = tfcomp2rotate_correctevent(TF_cal,TFs,transprop,spec_mat,NFFT,dt,npad0,npts,goodP,goodH,goodZ)

TF_name = TF_cal{1};

comp1 = TF_name(1);
comp2 = TF_name(4);
comp3 = TF_name(2);

TFidx12 = strmatch([comp1,comp2],{TFs.label},'exact');
TFidx32 = strmatch([comp3,comp2],{TFs.label},'exact');
TFidx13_2 = strmatch([comp1,comp3,'-',comp2],{TFs.label},'exact');

label_list{1} = [comp1,comp2];
label_list{2} = [comp1,comp3,'-',comp2];

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

if strcmp(comp3,'P')==1
    spec_3 = spec_mat(:,4);
    isgood3 = goodP;
else
    disp('Error, invalid component for rotation')
    return
end

hang = transprop.params.hang;
cang = cos(hang*pi/180);
sang = sin(hang*pi/180);

spec_2 = sang.*spec_H2+cang.*spec_H1;

corrspec1_2 = spec_1-((TFs(TFidx12).transfunc_tap)'.*spec_2);
corrspec3_2 = spec_3-((TFs(TFidx32).transfunc_tap)'.*spec_2);
corrspec1_32 = corrspec1_2-((TFs(TFidx13_2).transfunc_tap)'.*corrspec3_2);

amp1_2 = real(ifft(2*corrspec1_2,NFFT)./dt);
amp1_32 = real(ifft(2*corrspec1_32,NFFT)./dt);

corrtime1_2 = amp1_2(npad0+1:npad0+npts); 
corrtime1_32 = amp1_32(npad0+1:npad0+npts); 

corrgood1_2 = isgood1*isgood2;
corrgood1_32 = isgood1*isgood2*isgood3;




return