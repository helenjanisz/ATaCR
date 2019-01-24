function [lc2c1,label_list] = comp1rotate_calctransfunc(TF_cal,crr_stack)

% designed for either ZH or PH, removal of the horizontal component noise
% from either the Z or P component by rotating in the maximum direction of
% coherence

TF_name = TF_cal{1};
comp1 = TF_name(1);
comp2 = TF_name(2);

label_list{1}=[comp1,comp2];

if strcmp(comp1,'Z')==1
    gc2c1= crr_stack(:,3);
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

lc2c1 = gc2c1./gc2c2;

return