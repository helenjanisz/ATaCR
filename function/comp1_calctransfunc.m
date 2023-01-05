function [lc2c1,label_list,gc1c1_c2] = comp1_calctransfunc(TF_cal,cnn_stack,cnm_stack,f)
% function for calcating transfer function between two components
TF_name = TF_cal{1};
comp1 = TF_name(1);
comp2 = TF_name(2);

label_list{1}=[comp1,comp2]; 

if strcmp(comp2,'Z')==1
    gc2c2= cnn_stack(:,1);
    if strcmp(comp1,'1')==1
        gc2c1 = conj(cnm_stack(:,1));
    elseif strcmp(comp1,'2')==1
        gc2c1 = conj(cnm_stack(:,2));
    elseif strcmp(comp1,'P')==1
        gc2c1 = conj(cnm_stack(:,3));
    end
elseif strcmp(comp2,'1')==1
    gc2c2= cnn_stack(:,2);
    if strcmp(comp1,'Z')==1
        gc2c1 = cnm_stack(:,1);
    elseif strcmp(comp1,'2')==1
        gc2c1 = cnm_stack(:,4);
    elseif strcmp(comp1,'P')==1
        gc2c1 = cnm_stack(:,5);
    end
elseif strcmp(comp2,'2')==1
    gc2c2= cnn_stack(:,3);
    if strcmp(comp1,'Z')==1
        gc2c1 = cnm_stack(:,2);
    elseif strcmp(comp1,'1')==1
        gc2c1 = conj(cnm_stack(:,4));
    elseif strcmp(comp1,'P')==1
        gc2c1 = cnm_stack(:,6);
    end
elseif strcmp(comp2,'P')==1
    gc2c2= cnn_stack(:,4);
    if strcmp(comp1,'Z')==1
        gc2c1 = cnm_stack(:,3);
    elseif strcmp(comp1,'1')==1
        gc2c1 = conj(cnm_stack(:,5));
    elseif strcmp(comp1,'2')==1
        gc2c1 = conj(cnm_stack(:,6));
    end
end

if strcmp(comp1,'Z')==1
    gc1c1= cnn_stack(:,1);
elseif strcmp(comp1,'1')==1
    gc1c1= cnn_stack(:,2);
elseif strcmp(comp1,'2')==1
    gc1c1= cnn_stack(:,3);
elseif strcmp(comp1,'P')==1
    gc1c1= cnn_stack(:,4);
end 

lc2c1 = gc2c1./gc2c2;


cohc1c2 = abs(gc2c1).^2./(gc2c2.*gc1c1); %coherence between pressure and vertical
gc1c1_c2 = gc1c1.*(1-cohc1c2); %removing from the autospectra

% Debugging
% figure(99)
% clf
% loglog(f,gc1c1_c2);hold on
% loglog(f,gc1c1,'-r')

% figure(98)
% clf
% semilogx(f,cohc1c2_c2);hold on
% semilogx(f,cohc1c2,'-r');

return