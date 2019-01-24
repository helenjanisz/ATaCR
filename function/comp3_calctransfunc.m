function [lc2c1,lc2c3,lc2c4,lc3c4_c2,lc3c1_c2,lc4c1_c3c2,label_list] = comp3_calctransfunc(TF_cal,cnn_stack,cnm_stack,f)

TF_name = TF_cal{1};
comp1 = TF_name(1);
comp2 = TF_name(5);
comp3 = TF_name(4);
comp4 = TF_name(2);

label_list{1}=[comp1,comp2]; 
label_list{2}=[comp3,comp2];
label_list{3}=[comp4,comp2];
label_list{4}=[comp4,comp3,'-',comp2];
label_list{5}=[comp1,comp3,'-',comp2];
label_list{6}=[comp1,comp4,'-',comp3,comp2];

% define components
if strcmp(comp2,'Z')==1
    gc2c2= cnn_stack(:,1);
    if strcmp(comp1,'1')==1
        gc2c1 = conj(cnm_stack(:,1));
    elseif strcmp(comp1,'2')==1
        gc2c1 = conj(cnm_stack(:,2));
    elseif strcmp(comp1,'P')==1
        gc2c1 = conj(cnm_stack(:,3));
    end
    if strcmp(comp4,'1')==1
        gc2c4 = conj(cnm_stack(:,1));
    elseif strcmp(comp4,'2')==1
        gc2c4 = conj(cnm_stack(:,2));
    elseif strcmp(comp4,'P')==1
        gc2c4 = conj(cnm_stack(:,3));
    end
    if strcmp(comp3,'1')==1
        gc2c3 = conj(cnm_stack(:,1));
    elseif strcmp(comp3,'2')==1
        gc2c3 = conj(cnm_stack(:,2));
    elseif strcmp(comp3,'P')==1
        gc2c3 = conj(cnm_stack(:,3));
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
    if strcmp(comp4,'Z')==1
        gc2c4 = cnm_stack(:,1);
    elseif strcmp(comp4,'2')==1
        gc2c4 = cnm_stack(:,4);
    elseif strcmp(comp4,'P')==1
        gc2c4 = cnm_stack(:,5);
    end
    if strcmp(comp3,'Z')==1
        gc2c3 = cnm_stack(:,1);
    elseif strcmp(comp3,'2')==1
        gc2c3 = cnm_stack(:,4);
    elseif strcmp(comp3,'P')==1
        gc2c3 = cnm_stack(:,5);
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
    if strcmp(comp4,'Z')==1
        gc2c4 = cnm_stack(:,2);
    elseif strcmp(comp4,'1')==1
        gc2c4 = conj(cnm_stack(:,4));
    elseif strcmp(comp4,'P')==1
        gc2c4 = cnm_stack(:,6);
    end
    if strcmp(comp3,'Z')==1
        gc2c3 = cnm_stack(:,2);
    elseif strcmp(comp3,'1')==1
        gc2c3 = conj(cnm_stack(:,4));
    elseif strcmp(comp3,'P')==1
        gc2c3 = cnm_stack(:,6);
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
    if strcmp(comp4,'Z')==1
        gc2c4 = cnm_stack(:,3);
    elseif strcmp(comp4,'1')==1
        gc2c4 = conj(cnm_stack(:,5));
    elseif strcmp(comp4,'2')==1
        gc2c4 = conj(cnm_stack(:,6));
    end
    if strcmp(comp3,'Z')==1
        gc2c3 = cnm_stack(:,3);
    elseif strcmp(comp3,'1')==1
        gc2c3 = conj(cnm_stack(:,5));
    elseif strcmp(comp3,'2')==1
        gc2c3 = conj(cnm_stack(:,6));
    end
end

if strcmp(comp3,'Z')==1
    gc3c3= cnn_stack(:,1);
    if strcmp(comp1,'1')==1
        gc3c1 = conj(cnm_stack(:,1));
    elseif strcmp(comp1,'2')==1
        gc3c1 = conj(cnm_stack(:,2));
    elseif strcmp(comp1,'P')==1
        gc3c1 = conj(cnm_stack(:,3));
    end
    if strcmp(comp4,'1')==1
        gc3c4 = conj(cnm_stack(:,1));
    elseif strcmp(comp4,'2')==1
        gc3c4 = conj(cnm_stack(:,2));
    elseif strcmp(comp4,'P')==1
        gc3c4 = conj(cnm_stack(:,3));
    end
elseif strcmp(comp3,'1')==1
    gc3c3= cnn_stack(:,2);
    if strcmp(comp1,'Z')==1
        gc3c1 = cnm_stack(:,1);
    elseif strcmp(comp1,'2')==1
        gc3c1 = cnm_stack(:,4);
    elseif strcmp(comp1,'P')==1
        gc3c1 = cnm_stack(:,5);
    end
    if strcmp(comp4,'Z')==1
        gc3c4 = cnm_stack(:,1);
    elseif strcmp(comp4,'2')==1
        gc3c4 = cnm_stack(:,4);
    elseif strcmp(comp4,'P')==1
        gc3c4 = cnm_stack(:,5);
    end
elseif strcmp(comp3,'2')==1
    gc3c3= cnn_stack(:,3);
    if strcmp(comp1,'Z')==1
        gc3c1 = cnm_stack(:,2);
    elseif strcmp(comp1,'1')==1
        gc3c1 = conj(cnm_stack(:,4));
    elseif strcmp(comp1,'P')==1
        gc3c1 = cnm_stack(:,6);
    end
    if strcmp(comp4,'Z')==1
        gc3c4 = cnm_stack(:,2);
    elseif strcmp(comp4,'1')==1
        gc3c4 = conj(cnm_stack(:,4));
    elseif strcmp(comp4,'P')==1
        gc3c4 = cnm_stack(:,6);
    end
elseif strcmp(comp3,'P')==1
    gc3c3= cnn_stack(:,4);
    if strcmp(comp1,'Z')==1
        gc3c1 = cnm_stack(:,3);
    elseif strcmp(comp1,'1')==1
        gc3c1 = conj(cnm_stack(:,5));
    elseif strcmp(comp1,'2')==1
        gc3c1 = conj(cnm_stack(:,6));
    end
    if strcmp(comp4,'Z')==1
        gc3c4 = cnm_stack(:,3);
    elseif strcmp(comp4,'1')==1
        gc3c4 = conj(cnm_stack(:,5));
    elseif strcmp(comp4,'2')==1
        gc3c4 = conj(cnm_stack(:,6));
    end
end

if strcmp(comp4,'Z')==1
    gc4c4= cnn_stack(:,1);
    if strcmp(comp1,'1')==1
        gc4c1 = conj(cnm_stack(:,1));
    elseif strcmp(comp1,'2')==1
        gc4c1 = conj(cnm_stack(:,2));
    elseif strcmp(comp1,'P')==1
        gc4c1 = conj(cnm_stack(:,3));
    end
elseif strcmp(comp4,'1')==1
    gc4c4= cnn_stack(:,2);
    if strcmp(comp1,'Z')==1
        gc4c1 = cnm_stack(:,1);
    elseif strcmp(comp1,'2')==1
        gc4c1 = cnm_stack(:,4);
    elseif strcmp(comp1,'P')==1
        gc4c1 = cnm_stack(:,5);
    end
elseif strcmp(comp4,'2')==1
    gc4c4= cnn_stack(:,3);
    if strcmp(comp1,'Z')==1
        gc4c1 = cnm_stack(:,2);
    elseif strcmp(comp1,'1')==1
        gc4c1 = conj(cnm_stack(:,4));
    elseif strcmp(comp1,'P')==1
        gc4c1 = cnm_stack(:,6);
    end
elseif strcmp(comp4,'P')==1
    gc4c4= cnn_stack(:,4);
    if strcmp(comp1,'Z')==1
        gc4c1 = cnm_stack(:,3);
    elseif strcmp(comp1,'1')==1
        gc4c1 = conj(cnm_stack(:,5));
    elseif strcmp(comp1,'2')==1
        gc4c1 = conj(cnm_stack(:,6));
    end
end

lc2c1=gc2c1./gc2c2;
lc2c3=gc2c3./gc2c2;
lc2c4=gc2c4./gc2c2;

%coherences between same
cohc2c3=abs(gc2c3).^2./(gc2c2.*gc3c3);
cohc2c4=abs(gc2c4).^2./(gc2c2.*gc4c4);

%this is predicting the channels y, 2 and 3  from 1
gc3c3_c2=gc3c3.*(1-cohc2c3);
gc4c4_c2=gc4c4.*(1-cohc2c4);

%removing the effect of channel 1 from the cross spectra
gc3c1_c2=gc3c1-conj(lc2c3).*gc2c1;
gc4c1_c2=gc4c1-conj(lc2c4).*gc2c1;
gc3c4_c2=gc3c4-conj(lc2c3).*gc2c4;

%transfer function between 2 and 3 after removing the effect of channel 1
lc3c4_c2=gc3c4_c2./gc3c3_c2;
%transfer function between 2 and y after removing the effect of channel 1
lc3c1_c2=gc3c1_c2./gc3c3_c2;

%coherence between (2 and y), (3 and y) and (2 and 3) after removing effect
%of channel 1
cohc3c4_c2=abs(gc3c4_c2).^2./(gc4c4_c2.*gc3c3_c2);

gc4c4_c2c3=gc4c4_c2.*(1-cohc3c4_c2);
gc4c1_c2c3=gc4c1_c2-conj(lc3c4_c2).*gc3c1_c2;

lc4c1_c3c2 = gc4c1_c2c3./gc4c4_c2c3; % final TF

return

