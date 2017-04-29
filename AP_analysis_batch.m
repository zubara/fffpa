clear all;
close all;

%%CHANGE HERE
%%%%%
auto_file_list = 1; %1 for automatically addding all .*abf files in datadir
export = 'csv_batch'; % 'xls_batch' - 1 *.xls table per batch , 'xls'-1 *.xls file per cell, 'csv' - 1 csv file per cell

%%set path to files (workdir) and where to save outputs
datadir = '/u/62/zubarei1/unix/Desktop/elinas_APs/data_elina/';% workdir should ONLY contain analyzed files
savedir = '/u/62/zubarei1/unix/Desktop/elinas_APs/results/';%
%%%%
cd(datadir);

if auto_file_list == 0
    % type filenames here separated by commas here
    abflist = { '16411383.abf','16411379.abf',...	
                '16411457.abf',...	
%                 '16411143.abf','16411117.abf',...	
%                 '16411158.abf',  '16411181.abf',...
%                 '16411125.abf', '16411163.abf',...
%                 '16411054.abf',  '16411097.abf',...
%                 '2072.abf', '16411187.abf',...
%                 '2187.abf',...
%                 '16401015.abf', '16411088.abf',...
%                  '16411205.abf',...
%                 '16411120.abf'  
               };%
else
%%this is for automatically adding all files in the workdir into the
%%analysis (uncomment)
    workfiles = dir('*.abf');
    for n = 1:length(workfiles)
        abflist{n} = workfiles(n).name;
    end
end

% Prepare output
spikes = {length(abflist),1};
sweeps = cell(length(abflist),1);
cells = zeros(length(abflist),9);

cstats = cell(length(abflist),1);
X_2sw = zeros(length(abflist),29);
swstats = cell(length(abflist));

for f = 1:length(abflist)
    %%Modify this line accordingly
   [cstats{f},swstats{f}] = APanalysis(abflist{f},0.0697,0.896,'full',1,0);
   
   
   %%Arrange and save output
   sweeps_raw = zeros(length(swstats{f}),12);
   spikes_raw = zeros(length(swstats{f}),8);
           cells(f,:) = [cstats{f}.rmp;
                         cstats{f}.R_input;
                         cstats{f}.sag_dp; 
                         cstats{f}.rheobase_c;
                         cstats{f}.first_spike_sweep; 
                         cstats{f}.max_spike_sweep;
                         cstats{f}.rmp_std;
                         cstats{f}.tau0;
                         cstats{f}.tau1;]';
    
    for s = 1:length(swstats{f})
        if ~isempty(swstats{f}(s).summary)
            sweeps_raw(s,:) = [swstats{f}(s).count;
                               swstats{f}(s).current_step;  
                               swstats{f}(s).first_spike_lat;
                               swstats{f}(s).fmax_init;
                               swstats{f}(s).fmax_ss;
                               swstats{f}(s).adaptation_ratio;
                               swstats{f}(s).spont_spikes;
                               swstats{f}(s).post_rebound;
                               swstats{f}(s).adp_amp; 
                               swstats{f}(s).adp_lat;
                               swstats{f}(s).adp_theta;
                               swstats{f}(s).adp_norm]';
            spikes_raw(s,:) = mean(swstats{f}(s).summary,1);
        else
            sweeps_raw(s,:) = 0;
            spikes_raw(s,:) = 0;
        end
        spikes{f,s} = swstats{f}(s).summary;  
    end
    sweeps{f} = sweeps_raw;
    %% Pack 'compact' version for export to Matlab
    features_2sw = [cells(f,[1:4, 7:9])';
                    sweeps_raw(1,[1,3,9:12])'; %1st sweep with spikes
                    spikes_raw(1,[2,5:8])';%1st sweep with spikes
                    sweeps_raw(end,1:6)';%last sweep with max spikes
                    spikes_raw(end,[2,5:8])'];%last sweep with max spikes
    X_2sw(f,:) = features_2sw';
    
    cells_header = [{'Resting MP'}, {'Input_resistance'}, {'SAG DP'},...
                    {'Rheobase current'},{'1st Spike Sweep'},{'Max spike Sweep'},...
                    {'RMP std'},{'tau start'},{'tau stop'}];
                        
    sweeps_header = [{'Spike count'}, {'Current Step'},{'First Spike Lat'}...
                     {'Fmax Init'}, {'Fmax SS'}, {'Apadtation ratio'}, ...                            
                     {'Spont Spikes'}, {'Rebound'}, ...
                     {'ADP amplitude'},{'ADP latency'},{'ADP theta'},{'ADP norm'}];
                        
    spikes_header = [{'Peak Latency'}, {'AP amplitude'}, {'AP Rising Time'}, ...
                     {'AP Duration'}, {'AP HW'}, {'Thresh amplitude'},... 
                     {'AHP amplitude'}, {'AHP decay time'}];

    
    xls_header = [cells_header([1:4, 7:9])';
                  strcat('1_spk:',sweeps_header([1,3,9:12])'); %1st sweep with spikes
                  strcat('1_spk:',spikes_header([2,5:8])');%1st sweep with spikes
                  strcat('sat:',sweeps_header(1:6)');%last sweep with max spikes
                  strcat('sat:',spikes_header([2,5:8])')];%last sweep with max spikes

        %Write EXCEL file
if strcmp(export,'xls')
        xlswrite([savedir,abflist{f}(1:end-4),'.xlsx'], xls_header,'2Sweeps_data','A2');
        xlswrite([savedir,abflist{f}(1:end-4),'.xlsx'], abflist,'2Sweeps_data','B1')
        xlswrite([savedir,abflist{f}(1:end-4),'.xlsx'], features_2sw,'2Sweeps_data','B2')
elseif strcmp(export,'csv')
        fid = fopen([savedir,abflist{f}(1:end-4),'.csv'], 'w') ;
        for i = 1:length(xls_header)
            fprintf(fid, '%s,%d,\n', xls_header{i},features_2sw(i)) ;
        end
        fclose(fid) ;
end
 end

    if strcmp(export,'xls_batch')
            xlswrite([savedir,'abfdata',date,'.xlsx'], xls_header,'2Sweeps_data','A2');
            xlswrite([savedir,'abfdata',date,'.xlsx'], abflist,'2Sweeps_data','B1')
            xlswrite([savedir,'abfdata',date,'.xlsx'], X_2sw','2Sweeps_data','B2')
    elseif strcmp(export,'csv_batch')
%             fid = fopen([savedir,'abfdata',date,'.csv'], 'w') ;
%             for i = 1:length(xls_header)
%                 fprintf(fid, '%16s,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,\n', xls_header{i},X_2sw(:,i)) ;
%             end
%             fclose(fid) ;
            
    end


%%Export to MAT File
save([savedir,'abfdata_',date,'.mat'], 'abflist', 'spikes', 'sweeps',...
    'cells', 'X_2sw', 'cells_header', 'sweeps_header', 'spikes_header')
