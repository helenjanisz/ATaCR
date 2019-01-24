function [datout,new_taxis] = resample(datin,dt_new,dt_old);

if ~isempty(datin)
old_taxis = 0:dt_old:(length(datin)-1)*dt_old;
new_taxis = 0:dt_new:(length(datin)-1)*dt_old;
%apply anti-alias filter
fn = 1/2/dt_old;
fN = 1/2/dt_new;
w_c = fN/2;
w_c = w_c/fn;
[b,a] = butter(6,w_c,'low');
datin = filtfilt(b,a,datin);
datout = interp1(old_taxis,datin,new_taxis,'spline');
datout = datout';
new_taxis = new_taxis';
else
    datout = [];
    new_taxis = [];
end

return