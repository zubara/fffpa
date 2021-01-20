function [cell, traces]=APanalysis(filename, opt)
%Usage
%
%Input arguments:
%1) filename        - include the full name of the file with .abf in ' ... '
%2) current_start   - time when current step begins, in seconds from 
%                     the start of the sweep (e.g. 0.500 = 500 ms)
%3) current_stop    - time when current step ends, in seconds from the
%                     start of the sweep (e.g. 1.5 = 1500 ms)
%4) opt             - options structure

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
%   half_width     - AP Half-width: width halfway from thres to peak to the same
%                   voltage during repolarization (samples)
%   ahp_amp        - AfterHyperPolariaztion ap_thr - and the most negative 
%                   voltage during AHP (mV, size=number of spikes)  
%   ap_decay_time - Time interval between AP peak and AHP peak (in samples)
%   adp_amp        - Difference  between AHP peak and fast repolarization
%                   peak

% Written by Ivan Zubarev 2016 ivan.zubarev@aalto.fi

% Load data from .abf file. Using newer version abf2load:
[traces,si] = abfload(filename,'start',1,'stop','e');

%% *********************************************************************
cell = struct(); % parameters of the cell
cell.sweeps = struct();

current_start = opt.current_start;
current_stop = opt.current_stop;

cell.filename = filename; 
fs = 1e6/si; % Hz (check header)
t_scale = 1e-3*fs; %milliseconds
if t_scale ~= opt.t_scale
    disp('WARNING! Check sampling rate!!!')
    return
end
cell.t_scale = t_scale;
current_start = fix(current_start*fs);%fix(0.0697*fs);
cell.current_start = current_start;
current_stop  = fix(current_stop*fs);%fix(0.896*fs);


%Process opts
rb_window = opt.rb_window*t_scale;
not_a_train = opt.not_a_train*t_scale; % i.e. every 3 spikes that are more than 50ms away are considered separate.
volt_change_thres = opt.volt_change_thres; %i.e. threshold is detected as a first sample after which Voltage grows > 10 mV/ms
apr_thres = opt.apr_thres; %consider neighbouring sample within this threshold 'consequetive'
max_dep_rep = opt.max_dep_rep*t_scale; %maximium duration of depolarization/repolarization
last_adp = opt.last_adp*t_scale;

%%*********************************************************************

% if nargin < 5
%     spike_threshold = 1;
% end

% if nargin <6
%     viz = 1;
% end

totalSweeps = size(traces,3);   
spike_counts = zeros(totalSweeps,1);
spont_spikes = zeros(totalSweeps,1);
post_rebound = zeros(totalSweeps,1);



for sweep=1:totalSweeps    
        %%SPIKE DETECTION
        %Scan through all sweeps, count spikes, define 1st sweep with APs and 1st with most spikes
        [counts_old, indx_old, ~] = spike_times(traces(:,1,sweep), 1);

        [spike_counts(sweep), ~, spont_spikes(sweep), post_rebound(sweep)] = ...
            check_spont(counts_old, indx_old, current_start,current_stop,rb_window);   
end

%% TODO: apply indices, use soi
cell.first_rebound = find(post_rebound,1);
cell.first_AP = find(spike_counts,1); % 1st sweep with spikes

%[~, cell_stats.max_APS] = max(spike_counts); %1st sweep w maximum spikes
cell.max_APS = find(spike_counts==max(spike_counts),1,'last');


if  contains('all', opt.sweeps_flag) == 1
    %take all sweeps with spikes
    soif = cell.first_AP:cell.max_APS;
    sweeps_of_interest = find(spike_counts(soif))+cell.first_AP - 1;
    %disp(sweeps_of_interest)

else
    sweeps_of_interest = [];
    if contains('rhb', opt.sweeps_flag)
        sweeps_of_interest = [sweeps_of_interest, cell.first_AP];
    end
    if contains('sat', opt.sweeps_flag)
        sweeps_of_interest = [sweeps_of_interest, cell.max_APS];
    end 
        
end

cell.sweeps_of_interest = sweeps_of_interest;
cell.exported_sweeps = opt.default_export*ones(length(sweeps_of_interest),1);
if strcmp(opt.detect_adp, '2') == 1
    adp_k = length(sweeps_of_interest);
elseif strcmp(opt.detect_adp, '1') == 1
    adp_k = 1;
else
    adp_k = 0;
end
%disp(filename)
%disp(['Sweep # ',' #APs ', 'Rebound ', 'Spont '])
%disp([(1:totalSweeps)',spike_counts, post_rebound, spont_spikes])

%% Input resistance
%regression coefficient between the injected current and
%mean membrane voltage in sweeps starting from 100mV until first spike 
%computed over the last 50 ms of current step
if opt.inj_curr_channel == 1
    current_steps = round(squeeze(mean(traces(current_start:current_stop,2,:),1)));
elseif opt.inj_curr_channel == 0
    last_step = opt.inj_curr0+opt.inj_curr_incr*(size(traces,3) - 1);
    current_steps = (opt.inj_curr0:opt.inj_curr_incr:last_step)';
end
volt_mean = squeeze(mean(traces(current_stop-50*t_scale:current_stop,1,:))); %last 50 ms
[v_err,t_ind] = min(abs(volt_mean+100));
if v_err > mean(diff(current_steps))
  t_ind = 1;
end
if isempty(cell.first_AP)
  cell.first_AP = totalSweeps;    
end
if cell.first_AP == 1
  disp('Warning! APs detected in the first sweep => input resistance may be computed inaccurately!')
  b = glmfit(current_steps,volt_mean);
else
  b = glmfit(current_steps(t_ind:cell.first_AP-1),volt_mean(t_ind:cell.first_AP-1));
end
cell.R_input = b(2)*1e+3;
%% Rheobase

cell.rheobase_c = current_steps(cell.first_AP);%current step amplitude on which the 1st AP appear
%% Time constant
        
dt = traces(current_start:current_stop,1,t_ind);
t = (1./fs)*(0:(size(dt,1)-1))';
g = fittype('a-b*exp(-c*x)');
f0 = fit(t,dt,g,'StartPoint',[[ones(size(t)), -exp(-t)]\dt; 1]);
m0 = f0(t);
[~,tau_ind0] = min(abs((f0.a-m0)./f0.b-1/exp(1)));
cell.tau0 = tau_ind0/t_scale;

dt1 = traces(current_stop:end,1,t_ind);
t1 = (1./fs)*(0:(size(dt1,1)-1))';
f1 = fit(t1,dt1,g,'StartPoint',[[ones(size(t1)), -exp(-t1)]\dt1; 1]);
m1 = f1(t1);
[~,tau_ind1] = min(abs((f1.a-m1)./f1.b-1/exp(1)));
cell.tau1 = tau_ind1/t_scale;
        
%% Resting membrane potential (rmp)
%mean of MP's in the interval from 5 ms after sweep start to the 
%moment of current start across the sweeps
rmps  = squeeze(mean(traces(50:current_start,1,:),1));
cell.rmp = mean(rmps);
cell.rmp_std = std(rmps);
          
%% Sag depolarizeation
%ratio between maximum hyperpolaization in sweep with 100 mV 
%membrane voltage and mean voltage computed over last 50 ms of the
%current step
minhp = min(traces(:,1,t_ind));
av_ss = volt_mean(t_ind);
cell.sag_dp = minhp/av_ss;
disp(['SAG =', num2str(cell.sag_dp)]);

for j = 1:length(sweeps_of_interest)
    sweep = sweeps_of_interest(j);
    [count, peak_idx, peak_volt] = spike_times(traces(:,1,sweep),1);

% Now that the baseline stats are computed we get rid of everything outside the current step 
    if spont_spikes(sweep)==1 || post_rebound(sweep)==1
        peak_volt = peak_volt(logical((peak_idx>current_start).*(peak_idx<current_stop)));
        [count, peak_idx, ~,~] = check_spont(count, peak_idx, current_start,current_stop,rb_window);  
        traces([1:current_start-1,current_stop+max_dep_rep:end],1,sweep) = rmps(sweep);
    end

    if count > 0
%% Compute AP_thres, AP_amplitude, and AP_rising_time
      first_spike_lat = (peak_idx(1)-current_start);
      last_peak_events_tw = max((current_stop-peak_idx(end))/2,max_dep_rep);
      peak_isi = [diff(peak_idx);last_peak_events_tw];

      if peak_isi(1)/first_spike_lat > 5 && peak_isi(1) > median(peak_isi) && first_spike_lat < 30*t_scale
            %disp('DELAYED FIRING CELL WITH EARLY SPIKE!')
            if opt.reject_early_spikes == 1
                peak_idx = peak_idx(2:end);
                peak_isi = peak_isi(2:end);
                count = count-1;
                peak_volt = peak_volt(2:end);
            end
      end
          
      th_latency = zeros(count,1);
      ahp_v = zeros(count,1);
      ahp_ind = zeros(count,1);      

      [b,a] = butter(6,0.25,'low');
      cvd = filtfilt(b,a,traces(:,1,sweep));
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
              prev_thres = traces(th_latency(p-1),1,sweep);
              [~,tmp] = min(abs(traces(peak_idx(p):-1:peak_idx(p)-2*max_dep_rep,1,sweep)-prev_thres));
              tmp = peak_idx(p)-tmp;% - max_dep_rep;
              disp([' Threshold detection FAILED! Setting AP thres to a sample closest to threshold of previous spike! (spike #', num2str(p),')'])
          elseif isempty(tmp) && p == 1
              [~,tmp] = min(abs(traces(peak_idx(p):-1:peak_idx(p)-2*max_dep_rep,1,sweep)+20));
              tmp = peak_idx(p)-tmp;% - max_dep_rep;
              disp([' Threshold detection FAILED! Setting AP thres to a sample closest to -20uV! (spike #', num2str(p),')'])
          end
          %disp(tmp)
          dapind = find(-1*diff(tmp)>=apr_thres,1,'first'); %if == 1 means that the voltage is also growing fast in the next sample (for dealing with noise)
          if isempty(dapind)
              %disp(['DAPIND IS EMPTY! spike #', num2str(p)])
              %[~,dapind] = min(abs(data(tmp,1,sweep)+20)); % if (say) only one above threshold sample in AP rising, use closest ot 0mV
              dapind = length(tmp);
          end
          th_latency(p) = tmp(dapind);
          ap_indx(dep_ind) = [];
      end
      th_amp = traces(th_latency,1,sweep);
%         disp([size(th_latency),size(peak_volt)])
      if any(size(th_amp) ~= size(peak_volt))
          disp('OOPS! Something is wrong, this case needs attention.')
      end

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
          [ahp_v(i), ahp_ind(i)]= min(traces(peak_idx(i):ahp_rb(i),1,sweep));
      end
      
      ahp_ind = ahp_ind+peak_idx;


%% ADP computation
      %FIX definition of lask peak in the train
      %FIX adp window definition similarly to AHP
      %defines last peak in the first train (if ISI > not_a_train, or 
      %if isi is twice longer that the previous ISI)
      %finds respective AHP 
      %
      %adp_ind      -maximum depolarizatoin time (in samples)
      %adp_amp      -ADP amplitude measured from the respective AHP
      last_ahp = ahp_ind(end);
      if j <= adp_k % only computed for the first spike 
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
          if ahp_v(k)<=traces(next_event,1,sweep)
            adp2_bl = linspace(ahp_v(k),traces(next_event,1,sweep),next_event-last_ahp+1);
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

          else
              adp_wind = min(adp_wind(end),length(adp2_bl))-1;
              [adp_v,adp1_ind] = max(cvd(last_ahp:last_ahp+adp_wind));

            %  disp((adp_wind - adp2_ind)/(adp_wind - adp1_ind))
              if (adp_v-adp2_bl(adp1_ind)) < adp_rv - 1 || adp_rv < 1 || (adp_v-adp2_bl(adp1_ind))< 1
                  adp_v = 0;
                  adp_amp = 0;
                  adp_ind = last_ahp;                
              elseif (adp_wind - adp2_ind)/(adp_wind - adp1_ind+1) > 2.5

                  adp_v = 0;
                  adp_amp = 0;
                  adp_ind = last_ahp;                
              else
                  %%%CHOSE ADP
                adp_ind = adp2_ind+last_ahp;
                adp_v = cvd(adp_ind);
                adp_amp = adp_v-ahp_v(k);
              end
          end
      else
          adp_v = 0;
          adp_amp = 0;
          adp_ind = last_ahp; 
          %adp_lat = (adp_ind-last_ahp)./t_scale;
          %theta = adp_amp/(sqrt(adp_lat^2+adp_amp^2)+1e-6);
          %disp(['theta=', num2str(round(theta,2)),' norm=', num2str(sqrt(adp_lat^2+adp_amp^2))])

      end

      %% fmax_init, fmax_ss, and adaptation ratio
      %fmax_init    -inverse of shortest ISI between first 3 spikes
      %fmax_ss      -inverse of the mean ISI between last 4 spikes
      %adaptation_ratio - fmax_ss/fmax_init;   

      if  length(peak_isi) > 4
          fmax_init = max(1./(peak_isi(1:2)))*fs;
          fmax_ss = median((1./peak_isi(end-4:end-1)))*fs; %last peak_isi is the time from last peak to current end
      else 
          fmax_init = 1;
          fmax_ss = 1; %last peak_isi is the time from last peak to current end
          if sweep == cell.max_APS
            disp('Warning! could not compute setady state frequency (less than 3 spikes in the sweep!!)')
          end
      end
      adaptation_ratio = fmax_ss/fmax_init;        
    
        %% Sweep stats
        cell.sweeps(j).spont_spikes = spont_spikes(j);
        cell.sweeps(j).post_rebound = post_rebound(j);
        cell.sweeps(j).first_spike_lat = (peak_idx(1)-current_start); %./t_scale;
        cell.sweeps(j).current_step = current_steps(sweep);
        cell.sweeps(j).count = count;
        % if applicable
        cell.sweeps(j).fmax_ss = fmax_ss;
        cell.sweeps(j).fmax_init = fmax_init;
        cell.sweeps(j).adaptation_ratio = adaptation_ratio;

        %if applicable
        cell.sweeps(j).adp_lat = adp_ind-last_ahp;%./t_scale;
        cell.sweeps(j).adp_amp = adp_amp;
        cell.sweeps(j).last_ahp = last_ahp;%./t_scale;%for plotting only


        cell.sweeps(j).sweep = sweep;
        

        %% Spike stats
        %%plotting only
        cell.sweeps(j).ahp_lat = ahp_ind; %time instances of AHP  for plotting
        cell.sweeps(j).th_lat = th_latency; % threshold latency, samples
        cell.sweeps(j).peak_lat = peak_idx; % -current_start)./t_scale
        cell.sweeps(j).adp_ind = adp_ind;%./t_scale;


        %threshs
        cell.sweeps(j).th_amp = th_amp; % Voltage of detected threshold
        cell.sweeps(j).peak_volt = peak_volt;
        cell.sweeps(j).ahp_v = ahp_v;
        cell.sweeps(j).adp_v = adp_v;
            %peaks
    end
end
end
    %%
    
   



    function [count, indx, spont, rebound] = check_spont(counts_old, indx_old, current_start,current_stop,rb_window)    
    %Cheks if the cell is displaying spontaneous activity
            spont = 0;
            rebound = 0;
            if counts_old ~=0
                %check if all spikes are within current step
                sp_idxf = logical((current_start < indx_old) .* (indx_old<current_stop));
                if current_start < indx_old(find(sp_idxf==0,1)) < current_stop + rb_window
                    rebound = 1;
                end
                if nnz(sp_idxf==0) > 1
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
    
    
    function [N , idxs, volts] = spike_times(trace,threshold)
    above_thresh=find(trace(1:end) > threshold)  ;  
    if ~isempty(above_thresh)
        inddif = diff(above_thresh);
        rises = [above_thresh(1);above_thresh(find(inddif>3)+1)];
        falls = [above_thresh(inddif>3);above_thresh(end)];
        nspikes = length(rises);
        volts = zeros(nspikes,1);
        idxs = zeros(nspikes,1);
        for i = 1:nspikes
           [volts(i),idxs(i)] = max(trace(rises(i):falls(i)));
           idxs(i) = idxs(i) + rises(i)-1;
        end

    else
        volts=[];
        idxs = [];
        %disp('no spikes in trace')
    end

    N=length(volts) ;
    
    end
    
 
