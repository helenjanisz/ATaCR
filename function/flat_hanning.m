function win = flat_hanning(vec_t,tapertime)
% function to apply a box-car window with a hanning taper end.

win=ones(size(vec_t));

%ndata=length(vec_t);
dt=vec_t(2)-vec_t(1);


taperlength = floor(tapertime/dt);

taper = hanning(taperlength*2);
win(1:taperlength) = taper(1:taperlength);
win(end-taperlength+1:end) = taper(taperlength+1:end);

end
