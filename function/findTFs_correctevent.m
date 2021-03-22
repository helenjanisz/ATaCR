function trans_filename=correctevent_findTFs(inpath_trans,tf_op,station,eventid);

trans_mat_files = dir(fullfile(inpath_trans,sprintf('*.mat')));

if ~isempty(trans_mat_files)
    if tf_op == 2
        trans_filename = sprintf('%s%s_AVERAGE_transfun.mat',inpath_trans,station);
    elseif tf_op ==1
        eventid_num = datenum(eventid(1:8),'yyyymmdd');
        for it = 1:length(trans_mat_files)
            idx1 = strfind(trans_mat_files(it).name,'_');
            idx  = idx1(1);
            if ~isempty(idx)
                tfid(it,:) = trans_mat_files(it).name((idx+1):(idx+12));
                tfid_num(it) = datenum(trans_mat_files(it).name(idx+1:idx+8),'yyyymmdd');
            end
        end
        dday = eventid_num-tfid_num;
        iday = find(dday > 0); % find days prior to event
        [~,imin] = min(dday(iday)); % find min index
        min_idx = iday(imin); % index day nearest to event
        if isempty(min_idx) % if no available day prior, use nearest day
            [min_diff,min_idx]=min(abs(eventid_num-tfid_num)); %find tf from closest day to event
        end
        dayid = tfid(min_idx,:);
        trans_filename = sprintf('%s%s_%s_transfun.mat',inpath_trans,station,dayid);
    end
else
    trans_filename=nan;
end



return