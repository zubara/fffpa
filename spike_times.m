function [N , idxs, volts] = spike_times(trace,threshold1)

%   This function detects and locates the time points of action potentials in a trace of 
%   membrane potential as a function of time in a neuron. The trace should represent
%   a current clamp recording from a neuron.
%   Input: 
%   "trace" is the membrane voltage array of the neuron
%   "Theshold" is the value for the spike to cross to be detected.
%   Output:
%   The output array is the index location of spikes.
%
%   Rune W. Berg 2006
%   rune@berg-lab.net
%   www.berg-lab.net
%   Modified by Rune Berg May 2015

%   Modified by Ivan Zubarev 2016
%   ivan.zubarev@aalto.fi

set_crossgi=find(trace(1:end) > threshold1)  ;  % setting the threshold

if ~isempty(set_crossgi)
    inddif = diff(set_crossgi);
    set_cross_plusgi = [set_crossgi(1);set_crossgi(find(inddif>3)+1)];
    set_cross_minusgi = [set_crossgi(inddif>3);set_crossgi(end)];
    nspikes = length(set_cross_plusgi);
    volts = zeros(nspikes,1);
    idxs = zeros(nspikes,1);
    for i = 1:nspikes
       [volts(i),idxs(i)] = max(trace(set_cross_plusgi(i):set_cross_minusgi(i)));
       idxs(i) = idxs(i) + set_cross_plusgi(i)-1;
    end
    
else
    volts=[];
    idxs = [];
    %display('no spikes in trace')
end

N=length(volts) ;




