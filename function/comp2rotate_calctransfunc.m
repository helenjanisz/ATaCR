function [lc2c1,lc2c3,lc3c1_c2,label_list] = comp2rotate_calctransfunc(TF_cal,crr_stack,f)

TF_name = TF_cal{1};
comp1 = TF_name(1);
comp2 = TF_name(4);
comp3 = TF_name(2);

label_list{1}=[comp1,comp2]; 
label_list{2}=[comp3,comp2];
label_list{3}=[comp1,comp3,'-',comp2];

if strcmp(comp1,'Z')==1
    gc2c1= crr_stack(:,3);
    gc1c1= crr_stack(:,1);
else
    disp('Error, invalid component for rotation')
    return
end 

if strcmp(comp2,'H')==1
    gc2c2 = crr_stack(:,2);
else 
    disp('Error, invalid component for rotation')
    return
end

if strcmp(comp3,'P')==1
    gc2c3 = crr_stack(:,4);
    gc3c3 = crr_stack(:,5);
    gc3c1 = crr_stack(:,6);
else 
    disp('Error, invalid component for rotation')
    return
end

lc2c1 = gc2c1./gc2c2;
lc2c3 = gc2c3./gc2c2;

%coherences between same
cohc2c3=abs(gc2c3).^2./(gc2c2.*gc3c3);

%this is predicting the channels y, 2 and 3  from 1
gc3c3_c2=gc3c3.*(1-cohc2c3);

%removing the effect of channel 1 from the cross spectra
gc3c1_c2=gc3c1-conj(lc2c3).*gc2c1;

%transfer function between 2 and y after removing the effect of channel 1
lc3c1_c2=gc3c1_c2./gc3c3_c2;


% Debugging
% cohc1c2 = abs(gc2c1).^2./(gc2c2.*gc1c1); %coherence between pressure and vertical
% gc1c1_c2 = gc1c1.*(1-cohc1c2); %removing from the autospectra
% figure(99)
% clf
% loglog(f,gc1c1_c2);hold on
% loglog(f,gc1c1,'-r')
% pause

return