function plot_trace(ax, data, sweep)
    %ax = axes('Parent', f, 'Position', [0.13, 0.15, 0.775, 0.815]);%get(f, 'CurrentAxes');
    markers = {'bo', 'ro', 'co', 'ms'};
    tags = {'Peak', 'Thres', 'AHP', 'ADP'};
    %handles = {}
    ap_times = [sweep.peak_lat, sweep.th_lat, sweep.ahp_lat];
    ap_values = [sweep.peak_volt, sweep.th_amp, sweep.ahp_v];

    hold(ax, 'off')
    plot((1:length(data)),data,'PickableParts','all', 'HitTest','on', ...
        'tag', sprintf('data_trace'), 'Parent', ax);
    set(ax, 'YLim', [-100,120]);
    %set(ax, 'YLim', [-100,120], 'XLim',[50,length(data)-50]);
    hold(ax, 'on')
    for i = 1:size(ap_times,2)
        plot(ap_times(:,i), ap_values(:,i), markers{i}, 'PickableParts','all', 'HitTest','on', ...
            'Parent', ax, 'tag', sprintf(tags{i}));
    end
    if sweep.adp_lat ~= 0 && sweep.adp_ind ~= sweep.last_ahp
        plot(sweep.adp_ind, sweep.adp_v, markers{4}, 'PickableParts','all', 'HitTest','on', ...
            'Parent', ax, 'tag', sprintf(tags{4}));
    end
    
    legend(ax, [{'Trace'}, tags]);
    legend(ax, 'boxoff');
    drawnow
  end