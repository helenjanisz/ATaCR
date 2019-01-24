%tf_taper
%cosine taper for transfer function for specified high and low frequency
%limits
function TF_out = tf_taper(TF,f,hp,lp,taper_width);
taper_points = floor(taper_width*length(f));

if hp~=0
    idx_hp = max(find(f<=hp));
else
    idx_hp=1;
end
if lp~=0
    idx_lp=min(find(f>=lp));
else
    idx_lp=length(TF);
end

for ip=1:length(TF);
    if ip<idx_hp-taper_points
        TF_out(ip)=0;
    elseif ip>idx_lp+taper_points
        TF_out(ip)=0;
    elseif ip>=idx_hp && ip<=idx_lp
        TF_out(ip)=TF(ip);
    elseif ip>=idx_hp-taper_points && ip<idx_hp
        TF_out(ip)=TF(ip)*cos(((idx_hp-ip)/taper_points)*(pi/2));
    elseif ip<=idx_lp+taper_points && ip>idx_lp
        TF_out(ip)=TF(ip)*cos(((ip-idx_lp)/taper_points)*(pi/2));
    end
end

return