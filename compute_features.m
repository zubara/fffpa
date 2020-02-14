 function cell = compute_features(cell, traces, opt)  
    for j = 1: length(cell.sweeps)
        count = cell.sweeps(j).count;
        %% Prepare AHP amp, ap_amp
        cell.sweeps(j).ap_decay_time = cell.sweeps(j).ahp_lat - cell.sweeps(j).peak_lat;
        cell.sweeps(j).ahp_amp = cell.sweeps(j).ahp_v - cell.sweeps(j).th_amp;
        cell.sweeps(j).ap_amplitude = cell.sweeps(j).peak_volt - cell.sweeps(j).th_amp;
        cell.sweeps(j).ap_rising_time = cell.sweeps(j).peak_lat - cell.sweeps(j).th_lat;

        %% AP Half-width computation
        hw_volt = cell.sweeps(j).ap_amplitude/2 + cell.sweeps(j).th_amp;
        rep_hw_ind = zeros(count,1);
        dep_hw_ind = zeros(count,1);
        cell.sweeps(j).ap_duration = zeros(count,1);
     
        for i = 1:count
          peak_idx = cell.sweeps(j).peak_lat;
          depolarization = traces(peak_idx(i):-1:peak_idx(i) - ...
              opt.max_dep_rep*opt.t_scale,1,cell.sweeps(j).sweep);
          reporalization = traces(peak_idx(i):peak_idx(i) + ...
              opt.max_dep_rep*opt.t_scale,1,cell.sweeps(j).sweep);
          [~,rep_hw_ind(i)] = min(abs(reporalization - hw_volt(i)));
          [~,dep_hw_ind(i)] = min(abs(depolarization - hw_volt(i)));
          [~,md_ind] = min(abs(reporalization - cell.sweeps(j).th_amp(i)));
          md_ind = md_ind + peak_idx(i) - 1;
          cell.sweeps(j).ap_duration(i) = md_ind - cell.sweeps(j).th_lat(i);
        end

        cell.sweeps(j).ap_hw = rep_hw_ind + dep_hw_ind - 2; 
    end
end