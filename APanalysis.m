function [out_c,out]=APanalysis(filename,current_start,current_stop,soi,height,viz)
%Usage
%[cell_data,sweep_data]=APanalysis_v1_1(filename,current_start,current_stop,soi,height,viz)
%
%Input arguments:
%1) filename        - include the full name of the file with .abf in ' ... '
%2) current_start   - time when current step begins, in seconds from 
%                     the start of the sweep (e.g. 0.500 = 500 ms)
%3) current_stop    - time when current step ends, in seconds from the
%                     start of the sweep (e.g. 1.5 = 1500 ms)
%4) soi             - sweeps of interest, if empty then two sweeps are
%                     analyzed (1) - first sweep where spikes appear; 
%                     (2) - first sweep with most of the spikes; if 'full',
%                     all sweeps with spikes are analyzed ; if array of
%                     sweep indices then only the selected indices are 
%                     analyzed and plotted;
%5) height          - spike detection threshold in mV, +1 mV by default
%6) viz             - 1 to plot results for checking, disable when processing a large batch of files
%                     by default plots first sweep with spikes and last
%                     sweep with most spikes (if soi='full' or empty)
%Output
%Cell-specific output parameters (size = scalar)S
%   sweeps_of_interest - indices of sweeps_analyzed
%   R_input         - input resistence
%   sag_dp          - sag depolarization
%   rheobase_c      - rheobase current current step on which the 1st AP
%                     appears
%
%
%Sweep-specific features (size = length(soi))
%   sweep           - sweep index
%   count           - AP count
%   rmp             - resting state membrane potential 
%
%Spike-sepcific features (size = count)
%   peak_idx        - peak latencies (in samples)
%   peak_volt       - voltage at peaks (mV)
%
% 
%%        for all sweeps of interest compute outputs
%   th_amp          - AP threshold - defined as 1st time instance after which
%                   the voltage difference is growing greater than or equal 
%                   to 5 mV (mV, size=number of spikes in the sweep)
%   th_latency      - AP threshold latencies, (in samples, size=number 
%                   of spikes)
%   ap_amplitude    - AP amplitude, defined as difference between voltage 
%                   at peak and th_amp (mV, size=number of spikes)
%   ap_rising_time - AP rising time, time between reachin AP threshold 
%                   and the AP peak(ms??);
%   ap_duration    - time interval between ap_thr and the same voltage
%                   during repolarization (ms??)
%   half_width     - AP Half-width width halfway from thres to the same
%                   voltage during repolarization (samples)
%   ahp_amp        - AfterHyperPolariaztion ap_thr - and the most negative 
%                   voltage during AHP (mV, size=number of spikes)  
%   ahp_decay_time - Time interval between AP peak and AHP peak (in samples)
%   adp_amp        - Difference  between AHP peak and fast repolarization
%                   peak

% Written by Ivan Zubarev 2016 ivan.zubarev@aalto.fi

out = struct();
out_c = struct();
% Load data from .abf file. Using newer version abf2load:
[d,si] = abfload(filename,'start',1,'stop','e');
%seDefault input parameters
fs = 1e6/si; % Hz (check header)
t_scale = 1e-3*fs; %milliseconds
current_start = fix(current_start*fs);%fix(0.0697*fs);
current_stop  = fix(current_stop*fs);%fix(0.896*fs);
rb_window = 200*t_scale;
not_a_train = 75*t_scale; % i.e. every 3 spikes that are more than 50ms away are considered separate.
volt_change_thres = 10.; %i.e. threshold is detected as a first sample after which Voltage grows > 10 mV/ms
apr_thres = 5; %consider neighbouring sample within this threshold 'consequetive'
max_dep_rep = 4*t_scale; %maximium duration of depolarization/repolarization
last_adp = 75*t_scale;
if nargin < 5
    height = 1;
end

if nargin <6
    viz = 1;
end

totalSweeps = size(d,3);   
spike_counts = zeros(totalSweeps,1);
spont_spikes = zeros(totalSweeps,1);
post_rebound = zeros(totalSweeps,1);

for sweep=1:totalSweeps    
        %%SPIKE DETECTION
        %Scan through all sweeps, count spikes, define 1st sweep with APs and 1st with most spikes
        [counts_old, indx_old, ~] = spike_times(d(:,1,sweep),height);
%        spike_counts(sweep) = counts_old; testing purposes
        [spike_counts(sweep), ~, spont_spikes(sweep), post_rebound(sweep)] = ...
            check_spont(counts_old, indx_old, current_start,current_stop,rb_window);   
end
first_rebound = find(post_rebound,1);
first_spike = find(spike_counts,1); % 1st sweep with spikes
%[~, mx_sp_ind] = max(spike_counts); %1st sweep w maximum spikes
mx_sp_ind = find(spike_counts==max(spike_counts),1,'last');
if nargin < 4
%     if ~isempty(first_rebound) && first_rebound < first_spike
%         disp('Rebound Detected!')
%         sweeps_of_interest = [first_rebound,first_spike,mx_sp_ind];
%     else
        sweeps_of_interest = [first_spike,mx_sp_ind];
%    end
elseif strcmp(soi,'full') == 1
    %take all sweeps with spikes
    soif = first_spike:mx_sp_ind;
    sweeps_of_interest = find(spike_counts(soif))+first_spike - 1;
    %disp(sweeps_of_interest)
else
    sweeps_of_interest = soi;
end

%disp(filename)
%disp(['Sweep # ',' #APs ', 'Rebound ', 'Spont '])
%disp([(1:totalSweeps)',spike_counts, post_rebound, spont_spikes])

%% Input resistance
          %regression coefficient between the infected current steps and
          %mean membrane voltage computed over last 50 ms of current
          %injection in sweeps with no spikes + the first one with spikes
          
          current_steps = round(squeeze(mean(d(current_start:current_stop,2,:),1)));
          volt_mean = squeeze(mean(d(current_stop-500:current_stop,1,:)));
          [v_err,t_ind] = min(abs(volt_mean+100));
          if v_err > mean(diff(current_steps))
              t_ind = 1;
          end
          if isempty(first_spike)
              first_spike = totalSweeps;    
          end
          if first_spike == 1
              disp('Warning! APs detected in the first sweep => input resistance may be computed inaccurately!')
              b = glmfit(current_steps,volt_mean);
          else
              b = glmfit(current_steps(t_ind:first_spike-1),volt_mean(t_ind:first_spike-1));
          end
          input_resistance = b(2);
%% Time constant
        
        dt = d(current_start:current_stop,1,t_ind);
        t = (1./fs)*(0:(size(dt,1)-1))';
        g = fittype('a-b*exp(-c*x)');
        f0 = fit(t,dt,g,'StartPoint',[[ones(size(t)), -exp(-t)]\dt; 1]);
        m0 = f0(t);
        [~,tau_ind0] = min(abs((f0.a-m0)./f0.b-1/exp(1)));
        tau0 = tau_ind0/t_scale;
        
        dt1 = d(current_stop:end,1,t_ind);
        t1 = (1./fs)*(0:(size(dt1,1)-1))';
        f1 = fit(t1,dt1,g,'StartPoint',[[ones(size(t1)), -exp(-t1)]\dt1; 1]);
        m1 = f1(t1);
        [~,tau_ind1] = min(abs((f1.a-m1)./f1.b-1/exp(1)));
        tau1 = tau_ind1/t_scale;
        
%% Resting membrane potential (rmp)
          %mean of MP's in the interval from 5 ms after sweep start to the 
          %moment of current start across the sweeps
          rmps  =squeeze(mean(d(50:current_start,1,:),1));
          rmp = mean(rmps);
          rmp_std = std(rmps);
          
%% Sag depolarizeation
          %ratio between maximum hyperpolaization in the very first sweep
          %and mean membrane voltage computed over last 50 ms of current
          %injection across all sweeps
          minhp = min(d(:,1,1));
          av_ss = volt_mean(1);
          sagdp = minhp/av_ss;
          disp(['SAG =', num2str(sagdp)]);

for j = 1:length(sweeps_of_interest)
    sweep = sweeps_of_interest(j);
    [count, peak_idx, peak_volt] = spike_times(d(:,1,sweep),height);

% Now that the baseline stats are computed we get rid of everything outside the current step 
    if spont_spikes(sweep)==1 || post_rebound(sweep)==1
        peak_volt = peak_volt(logical((peak_idx>current_start).*(peak_idx<current_stop)));
        [count, peak_idx, ~,~] = check_spont(count, peak_idx, current_start,current_stop,rb_window);  
        d([1:current_start-1,current_stop+max_dep_rep:end],1,sweep) = rmps(sweep);
    end

    if count > 0
%% Compute AP_thres, AP_amplitude, and AP_rising_time
          first_spike_lat = (peak_idx(1)-current_start);
          last_peak_events_tw = max((current_stop-peak_idx(end))/2,max_dep_rep);
          peak_isi = [diff(peak_idx);last_peak_events_tw];
                   
%           disp(first_spike_lat)
%           disp(peak_isi(1))
%           disp(median(peak_isi))
          if peak_isi(1)/first_spike_lat > 5 && peak_isi(1) > median(peak_isi)
                %disp('DELAYED FIRING CELL WITH EARLY SPIKE!')
                peak_idx = peak_idx(2:end);
                peak_isi = peak_isi(2:end);
                count = count-1;
                peak_volt = peak_volt(2:end);
          end
          
          th_latency = zeros(count,1);
          ap_durations = zeros(count,1);
          rep_hw_ind = zeros(count,1);
          dep_hw_ind = zeros(count,1);
          ahp_v = zeros(count,1);
          ahp_ind = zeros(count,1);      
        
          [b,a] = butter(6,0.25,'low');
          cvd = filtfilt(b,a,d(:,1,sweep));
          volt_diff=diff(cvd); %voltage difference with next time sample
          
          b = ones(t_scale,1); %sliding sum filter
          vd = conv(volt_diff,b,'same'); %low-pass-filtered version of the signal derivative used for detection of low frequency events
          
          ap_indx = find(vd>=volt_change_thres); %find indices of all samples where the voltage grows above the threshold
          for p =1:(count)
              dep_ind = ap_indx<peak_idx(p); %for ap_indx take only samples before the current p'th peak
              tmp = ap_indx(dep_ind);
              tmp = tmp(end:-1:1);%%tmp should be ordered vice versa for better threshold
              %%detection
              if isempty(tmp) || isempty(dep_ind)
                  disp(['WARNING! Cannot detect voltage change above threshold, trying lower threshold for spike #',num2str(p)])
                  ap_indx1 = find(vd(peak_idx(p-1):1:peak_idx(p))>=volt_change_thres*0.5);%try peak_idx(p)-max_dep_rep
                  %dep_ind1 = ap_indx1<peak_idx(p); %for ap_indx take only samples before the current p'th peak
                  tmp = ap_indx1+peak_idx(p-1);
                  tmp = tmp(end:-1:1);%%tmp should be ordered vice versa for better threshold
              end
              if isempty(tmp) && p > 1 % handle situations if nothing is above threshold 
                  prev_thres = d(th_latency(p-1),1,sweep);
                  [~,tmp] = min(abs(d(peak_idx(p):-1:peak_idx(p)-2*max_dep_rep,1,sweep)-prev_thres));
                  tmp = peak_idx(p)-tmp;% - max_dep_rep;
                  disp([' Threshold detection FAILED! Setting AP thres to a sample closest to threshold of previous spike! (spike #', num2str(p),')'])
              elseif isempty(tmp) && p == 1
                  [~,tmp] = min(abs(d(peak_idx(p):-1:peak_idx(p)-2*max_dep_rep,1,sweep)+20));
                  tmp = peak_idx(p)-tmp;% - max_dep_rep;
                  disp([' Threshold detection FAILED! Setting AP thres to a sample closest to -20uV! (spike #', num2str(p),')'])
              end
              %disp(tmp)
              dapind = find(-1*diff(tmp)>=apr_thres,1,'first'); %if == 1 means that the voltage is also growing fast in the next sample (for dealing with noise)
              if isempty(dapind)
                  %disp(['DAPIND IS EMPTY! spike #', num2str(p)])
                  %[~,dapind] = min(abs(d(tmp,1,sweep)+20)); % if (say) only one above threshold sample in AP rising, use closest ot 0mV
                  dapind = length(tmp);
              end
              th_latency(p) = tmp(dapind);
              ap_indx(dep_ind) = [];
          end
          th_amp = d(th_latency,1,sweep);
 %         disp([size(th_latency),size(peak_volt)])
          if any(size(th_amp) ~= size(peak_volt))
              disp('OOPS! Something is wrong, this case needs attention.')
          end
          
          ap_amplitude = peak_volt - th_amp;
          ap_rising_time = peak_idx - th_latency;
          
%% AP Half-width computation
          hw_volt = ap_amplitude/2+th_amp;
          for i = 1:count
              depolarization = d(peak_idx(i):-1:peak_idx(i)-max_dep_rep,1,sweep);
              reporalization = d(peak_idx(i):peak_idx(i) + max_dep_rep,1,sweep);
              [~,rep_hw_ind(i)]=min(abs(reporalization - hw_volt(i)));
              [~,dep_hw_ind(i)]=min(abs(depolarization - hw_volt(i)));
              [~,md_ind]=min(abs(reporalization - th_amp(i)));
              md_ind = md_ind+peak_idx(i)-1;
              ap_durations(i) = md_ind - th_latency(i);
          end
          half_width = rep_hw_ind+dep_hw_ind-2; 
%% AHP computation   
          ahp_rb = zeros(count,1);
         
          for i = 1:count %minimum between peaks
              % find first 10 time instances where lp-filtered signal goes up
              % after the peak
              ap_rb_ind1 = find(volt_diff(peak_idx(i):fix(peak_idx(i)+peak_isi(i)/2))>0,10);  
              if isempty(ap_rb_ind1) %|| ap_rb_ind1(end)<max_dep_rep
                  %disp('Searching AHPS in entire ISI/2')
                  %if it doesn't rais warning and look for minimum at 1/2
                  %of ISI
                  ahp_rb(i) = fix(peak_idx(i)+peak_isi(i)/2);
                 % disp(ahp_rb(i)-peak_idx(i))
              else
                  %take twice as much time window and look for minimum
                  %between it an the peak
                  ap_rb_ind1 =  ap_rb_ind1(end);
                  ahp_rb(i) = fix(peak_idx(i)+min(ap_rb_ind1*2,fix(peak_isi(i)/2)));
              end
              
              [ahp_v(i), ahp_ind(i)]= min(d(peak_idx(i):ahp_rb(i),1,sweep));
          end
          ahp_decay_time = ahp_ind;
          ahp_ind = ahp_ind+peak_idx;
          ahp_amp = ahp_v-th_amp;
          
%% ADP computation
          %FIX definition of lask peak in the train
          %FIX adp window definition similarly to AHP
          %defines last peak in the first train (if ISI > not_a_train, or 
          %if isi is twice longer that the previous ISI)
          %finds respective AHP 
          %
          %adp_ind      -maximum depolarizatoin time (in samples)
          %adp_amp      -ADP amplitude measured from the respective AHP
          if j == 1
              for k = 1:length(peak_isi)
                  if k < 2
                      if peak_isi(k)>=not_a_train
                         %last_peak = peak_idx(k);
                         break
                      end
                  else
                      if peak_isi(k)>=not_a_train || peak_isi(k)>=peak_isi(k-1)*3
                        break
                      end
                  end
              end
              last_ahp = ahp_ind(k);
              if k < count
                  next_event = th_latency(k+1);
                  adp_window = last_ahp + fix((next_event - last_ahp)/2);
                  
              else
                  adp_window = last_ahp+last_adp;
                  next_event = last_ahp+fix(peak_isi(k));
              end
              
              %find 10 instances where the lp-filtered signal goes down
              adp_wind = find(vd(last_ahp:adp_window)<-0.75,10,'last');               
              %find the maximum paek relative to the line from last ahp to
              %next event
              if ahp_v(k)<=d(next_event,1,sweep)
                adp2_bl = linspace(ahp_v(k),d(next_event,1,sweep),next_event-last_ahp+1);
              else
                  adp2_bl = repmat(ahp_v(k),next_event-last_ahp+1,1)';
              end
              
              
              [adp_rv,adp2_ind] = max(cvd(last_ahp:next_event) - adp2_bl');
              
              %disp(min(vd(last_ahp:adp_window)))
              %disp(std(vd(last_ahp:adp_window)))
              if isempty(adp_wind)
                  disp('ADP wind is empty!')
                  %adp_wind = adp_window - last_ahp;
                  adp_v = 0;
                  adp_amp = 0;
                  adp_ind = last_ahp;
                  adp1_ind = 0;
                  adp_wind = 0;
              else
                  adp_wind = min(adp_wind(end),length(adp2_bl))-1;
                  [adp_v,adp1_ind] = max(cvd(last_ahp:last_ahp+adp_wind));
                  if (adp_v-adp2_bl(adp1_ind)) < adp_rv - 1 || adp_rv < 1 || (adp_v-adp2_bl(adp1_ind))< 1
                      adp_v = 0;
                      adp_amp = 0;
                      adp_ind = last_ahp;                
                  else
                    adp_ind = adp1_ind+last_ahp;
                    adp_amp = adp_v-ahp_v(k);
                  end
                  disp(['ADP1=', num2str(round(adp_v-adp2_bl(adp1_ind),2)),' ADP2=', num2str(round(adp_rv,2))])
                  disp(['ADP1 lat=', num2str(round(adp1_ind),2),' ADP2 lat=', num2str(round(adp2_ind,2))])
              end
                %plot(adp2_ind,d(adp2_ind,1,sweep),'kx','markerSize',10)
              
              figure;hold on;plot((1:size(d,1)),d(:,1,sweep));
              plot((last_ahp:last_ahp+adp_wind),cvd(last_ahp:last_ahp+adp_wind),'r')
              plot(last_ahp:next_event,adp2_bl,'g--');
              if adp_v ~= 0
                  plot(last_ahp+adp2_ind,cvd(last_ahp+adp2_ind),'md')
                  plot(last_ahp+adp1_ind,d(last_ahp+adp1_ind,1,sweep),'go','MarkerSize',11)
              end
              
              disp('!')
              
              adp_lat = (adp_ind-last_ahp)./t_scale;
              %disp(['ADP amp=',num2str(adp_rv),' lat=',num2str(adp_lat)]);
              theta = adp_amp/(sqrt(adp_lat^2+adp_amp^2)+1e-6);
              %disp(['theta=', num2str(round(theta,2)),' norm=', num2str(sqrt(adp_lat^2+adp_amp^2))])
              
          end
          
%          disp('!')
          %[adp_v,adp_ind] = max(cvd(last_ahp:adp_window));
          

          
%% fmax_init, fmax_ss, and adaptation ratio
          %fmax_init    -inverse of shortest ISI between first 3 spikes
          %fmax_ss      -inverse of the mean ISI between last 4 spikes
          %adaptation_ratio - fmax_ss/fmax_init;   
          
          if  length(peak_isi) > 4
              fmax_init = max(1./(peak_isi(1:2)))*fs;
              fmax_ss = median((1./peak_isi(end-4:end-1)))*fs; %last peak_isi is the time from last peak to current end
              adaptation_ratio =fmax_ss/fmax_init;        
          else
              fmax_init = NaN;
              fmax_ss = NaN;
              adaptation_ratio = NaN;
              if sweep == mx_sp_ind;
                disp('Warning! could not compute setady state frequency (less than 3 spikes in the sweep!!)')
              end
          end


    %end
%disp(current_steps);

%% Cell stats
        out_c.rmp = rmp;
        out_c.R_input  = input_resistance;
        out_c.sag_dp = sagdp; 
        out_c.rheobase_c   = current_steps(first_spike);%current step amplitude on which the 1st AP appear
        out_c.first_spike_sweep = first_spike;
        out_c.max_spike_sweep = mx_sp_ind;
        out_c.rmp_std = rmp_std;
        out_c.first_rebound = first_rebound;
        out_c.tau0 = tau0;
        out_c.tau1 = tau1;
        %out_c.last_rebound = last_rebound;
        %% Sweep stats
        out(j).spont_spikes = spont_spikes(j);
        out(j).post_rebound = post_rebound(j);
        out(j).first_spike_lat = (peak_idx(1)-current_start)./t_scale;
        out(j).current_step = current_steps(sweep);
        out(j).count = count;
        % if applicable
        out(j).fmax_ss = fmax_ss;
        out(j).fmax_init = fmax_init;
        out(j).adaptation_ratio = adaptation_ratio;
        %if applicable
        out(j).adp_lat = (adp_ind-last_ahp)./t_scale;
        out(j).adp_amp = adp_amp;
        out(j).adp_theta = theta;
        out(j).adp_norm = sqrt(adp_lat^2+adp_amp^2);
        out(j).last_ahp = last_ahp./t_scale;%for plotting only

        out(j).adp_v = adp_v;
        out(j).sweep = sweep;

        %% Spike stats
        %%plotting only
        out(j).ahp_lat = ahp_ind/t_scale; %time instances of AHP  for plotting
        out(j).ahp_rb = ahp_rb; % time window where to look for AHPs
        out(j).th_lat = (th_latency-current_start)./t_scale;
        %ahps
        out(j).ahp_amp = ahp_amp; % voltage difference at AHP vs AP threshold
        out(j).ahp_decay_time = ahp_decay_time./t_scale; % AHP latency from the corresponding spike 
        %threshs
        out(j).th_amp = th_amp; % Voltage of detected threshold
        %peaks
        out(j).peak_lat = (peak_idx-current_start)./t_scale;
        out(j).ap_amplitude = ap_amplitude;
        out(j).ap_rising_time = ap_rising_time./t_scale;
        out(j).ap_duration = ap_durations./t_scale;
        out(j).ap_hw = half_width./t_scale;
        out(j).header = [{'Peak Latency'}, {'AP amplitude'}, {'AP Rising Time'}, {'AP Duration'},...
                          {'AP HW'}, {'Thresh amplitude'}, {'AHP amplitude'},...
                          {'AHP decay time'}];
        %AP rising time variance? Duratio ~2 times RT, HW~ 0.5 of Duration, Th AMP
        %growing AHP ampp growing, AHP decay time saturatings
        % duration and growing?
        out(j).summary = [out(j).peak_lat,out(j).ap_amplitude,out(j).ap_rising_time, ...
                            out(j).ap_duration, out(j).ap_hw,out(j).th_amp, ...
                            out(j).ahp_amp,out(j).ahp_decay_time];
        %     else
        %         out(j)
    end 

end
if viz == 1 && ~isempty(sweeps_of_interest)
    cs_sc = current_start./t_scale;
    f=figure('units','normalized','outerposition',[0 0 1 1]);
    %% FIX plotting for 'full'
    suptitle(filename)
    if nargin < 2 || strcmp(soi,'full') == 1 || length(soi)==2;
        subplot(1,2,1)
        sweep = sweeps_of_interest(1);
        plot((1:size(d,1))./t_scale,d(:,1,sweep));hold on; 
        plot(out(1).peak_lat+cs_sc,d(fix(out(1).peak_lat*t_scale)+current_start,1,sweep),'bo');
        plot(out(1).peak_lat-out(1).ap_rising_time+cs_sc,out(1).th_amp,'ro');%hold off;
        plot(out(1).peak_lat+out(1).ahp_decay_time+cs_sc,out(1).ahp_amp+out(1).th_amp,'co');
        if adp_amp~=0
            plot(out(1).adp_lat+out(1).last_ahp,out(1).adp_v, 'ms');
        end
        title(['First sweep with spikes #', num2str(sweep)]);
        legend('Data', 'AP peaks', 'AP thres', 'AHP', 'ADP');

        subplot(1,2,2)
        sweep = sweeps_of_interest(end);
        plot((1:size(d,1))./t_scale,d(:,1,sweep));hold on; 
        plot(out(end).peak_lat+cs_sc,d(fix(out(end).peak_lat*t_scale)+current_start,1,sweep),'bo');
        plot(out(end).peak_lat-out(end).ap_rising_time+cs_sc,out(end).th_amp,'ro');%hold off;
        plot(out(end).peak_lat+out(end).ahp_decay_time+cs_sc,out(end).ahp_amp+out(end).th_amp,'co');
        legend('Data', 'AP peaks', 'AP thres', 'AHP');
        title(['First sweep with most spikes #',num2str(sweep)]);
        hold off;
     else
        dispitems = min(6,length(sweeps_of_interest));
        for s =1:dispitems
            subplot(dispitems,1,s);
            sweep = sweeps_of_interest(s);
            plot((1:size(d,1))./t_scale,d(:,1,sweep));hold on; 
            plot(out(s).peak_lat+cs_sc,d(fix(out(s).peak_lat*t_scale)+current_start,1,sweep),'bo');
            plot(out(s).peak_lat-out(s).ap_rising_time+cs_sc,out(s).th_amp,'ro');
            plot(out(s).peak_lat+out(s).ahp_decay_time+cs_sc,out(s).ahp_amp+out(s).th_amp,'co');
            %plot(out(s).adp_lat,out(s).adp_v, 'ms');
            title(['Selected sweep #', num2str(sweep)]);
            %legend('Data', 'AP peaks', 'AP thres', 'AHP');
            %hold off;
        end
    end
end
%saveas(f,[filename(1:end-4),'.jpg'])
end


    function [count, indx, spont, rebound] = check_spont(counts_old, indx_old, current_start,current_stop,rb_window)    
            spont = 0;
            rebound = 0;
            if counts_old ~=0
                %check if all spikes are within current step
                sp_idxf = logical((current_start < indx_old) .* (indx_old<current_stop));
                if current_start < indx_old(find(sp_idxf==0,1)) < current_stop + rb_window
                    rebound = 1;
                end
                if nnz(sp_idxf==0) > 1;
                    spont = 1;
                end       
            end
            
            %leave only spikes within current step and rebounds
            if spont == 1 || rebound == 1
                indx = indx_old(sp_idxf);
                count = nnz(sp_idxf);
            else
                indx = indx_old;
                count = counts_old;
            end
    end
    

