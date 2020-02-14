function [header, results] = prepare_sweep_results(cell_obj, sweep_header, spike_header, opt)

    % works on a single cell
    % adjusts time scale for output
    % integrates over single spike data within a sweep (mean, 1st, median, variance)
    % saves stuff in a specified format assigning header names
    %disp(length(cell_obj.sweeps))
    %spike_fun = opt.spike_fun;
    %while size(spike_fun,1) < length(cell_obj.sweeps)
    %    spike_fun = [spike_fun;opt.spike_fun(end,:)];
    %end
    spike_fun = [opt.spike_fun(1), repmat(opt.spike_fun(2), 1, length(cell_obj.sweeps)-1)];
    time_scaled_vars = {'first_spike_lat', 'adp_lat', 'ap_decay_time', ...
                        'ap_rising_time', 'ap_duration', 'ap_hw'};

    %disp(size(spike_fun))
    
    exports = logical(cell_obj.exported_sweeps);
    if contains('all',opt.sweeps_flag)
        sweep_names = cell(sum(exports), 1);
        for zz = 1:sum(exports)
            sweep_names{zz} = ['sw', num2str(zz)];
        end
    else
        sweep_names = opt.sweeps_flag; 
    end
    sweep_data = zeros(sum(exports),length(sweep_header));
    spike_stats = zeros(sum(exports),length(spike_header));
    
    z = 1;    
    for i = 1:length(cell_obj.sweeps)
        if exports(i) == 1
            for jj = 1:length(sweep_header)
                if ismember(sweep_header{jj},time_scaled_vars) == 1
                    sweep_data(z,jj) = cell_obj.sweeps(i).(sweep_header{jj})./opt.t_scale;
                else
                    sweep_data(z,jj) = cell_obj.sweeps(i).(sweep_header{jj});    
                end
            end

            spike_data = zeros(length(spike_header), cell_obj.sweeps(i).count);
            for h = 1:length(spike_header)
                if ismember(spike_header{h},time_scaled_vars) == 1
                    spike_data(h,:) = cell_obj.sweeps(i).(spike_header{h})./opt.t_scale;
                else
                    spike_data(h,:) = cell_obj.sweeps(i).(spike_header{h});
                end
            end

            switch spike_fun(i)
                case '1' % take first
                    %disp('Taking First')
                    spike_stats(z,:) = spike_data(:,1)';
                case '2' % mean
                    %disp('Taking Mean')
                    spike_stats(z,:) = mean(spike_data,2)';
                case '3' % median
                    spike_stats(z,:) = median(spike_data,2)';
                case '4' %variance
                    spike_stats(z,:) = var(spike_data,2)';
                case '5' % mean diff
                    spike_stats(z,:) = mean(diff(spike_data,1,2),2)';
            end
            z = z + 1;
        end
    end
    results = [sweep_data, spike_stats]; % n_exports x n_features array
    %results = results(exports,:);
    header = [];
    new_header = [sweep_header, spike_header];
    %disp(length(cell_obj.sweeps))
    for zz = 1:sum(exports)
        nh = cell(1, length(new_header));%{'',''}
        for k = 1:length(new_header)
            nh{k} = [sweep_names{zz},'_', new_header{k}];
        end
        header = [header,nh];
    end
        %header = [header; new_header];
    results = results';
    results = results(:);

%     disp(size(results))
%     disp(size(header))
end