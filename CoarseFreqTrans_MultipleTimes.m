% This version of code is used to generate Symbol Error plots along
% desiredr frwquency band. For every frequency offset, Symbol error is
% calculated trial times and average result is recorded.

%   Özün Yiğit Bayram 260206036
%   Yahya Ekin 260206054


fofset = -65:5:65; trial = 10;
Symbol_Error_per_offset = zeros(1,length(fofset));
Freq_Estimation_per_offset = zeros(1,length(fofset));

Symbol_Error_Trial = zeros(1,trial);
for j=1:1:length(fofset)
freq_offset = fofset(j);
% pulrecsig.m: create pulse shaped received signal

for k=1:1:trial

N=100; M=20; Ts=.0001;   % # symbols, oversampling factor
time=Ts*N*M; t=Ts:Ts:time; % sampling interval & time vector
message=pam(N,4,5);              % 4-level signal of length N
mup=zeros(1,N*M); 
mup(1:M:N*M)=message;            % oversample by integer length M
ps=hamming(M);             % blip pulse of width M
s=filter(ps,1,mup);        % convolve pulse shape with data
fc=1000+freq_offset; phoff=0;       % carrier freq. and phase
c=cos(2*pi*fc*t+phoff);    % construct carrier
rsc=s.*c;                  % modulated signal (small carrier)

% Adding Noise Due to AWGN Channel
noisy_rsc = awgn(rsc,30,"measured");

% Receiver Side
r=noisy_rsc;                           % suppressed carrier r
q=r.^2;                                % square nonlinearity
fl=500; ff=[0 .38 .39 .41 .42 1];      % BPF center frequency at .4
fa=[0 0 1 1 0 0];                      % which is twice f_0
h=firpm(fl,ff,fa);                     % BPF design via firpm
rp=filter(h,1,q);                      % filter gives preprocessed r

% pllpreprocess.m: recover "unknown" freq and phase using FFT
fftrBPF=fft(rp);                       % spectrum of rBPF
[m,imax]=max(abs(fftrBPF(1:end/2)));   % find freq of max peak
ssf=(0:length(rp))/(Ts*length(rp));    % frequency vector
freqS=ssf(imax);                       % freq at the peak
phasep=angle(fftrBPF(imax));           % phase at the peak
[IR,f]=freqz(h,1,length(rp),1/Ts);     % freq response of filter
[mi,im]=min(abs(f-freqS));             % freq where peak occurs
phaseBPF=angle(IR(im));                % < of BPF at peak freq
phaseS=mod(phasep-phaseBPF,pi);        % estimated angle

% Demodulation of the signal
c2=cos(2*pi*(freqS/2)*t+phaseS/2);
demod_signal = noisy_rsc.*c2;          % demod received signal
Freq_Estimation_per_offset(j)=freqS/2;

% Low Pass Filter for Demodulated Signal
fc1 = freqS/2; 
[b1,a1] = butter(14,(fc1/((1/Ts)/2)),"low");
LWPfilteredSignal = 2*filter(b1,a1,demod_signal);

% Compensate for the delay
groupDelay = mean(grpdelay(b1, a1));   % Calculate group delay
LWPfilteredSignal = LWPfilteredSignal(groupDelay+1:end);

% Adjust the time vector for plotting
t_adjusted = t(1:length(LWPfilteredSignal));

% Cross-Correlation and Symbol Decision
y=xcorr(LWPfilteredSignal,ps);         % correlate pulse with received signal
z=y(N*M:M:2*N*M-1)/(pow(ps)*M);        % downsample to symbol rate and normalize
mprime=quantalph(z,[-3,-1,1,3])';      % quantize to +/-1 and +/-3 alphabet

% Error Calculation
error=0;
for i=1:length(message)
    if message(i) ~= mprime(i)
        error=error+1;
    end
end
%disp("Normalized Symbol Error:"); 
p1 = (error/length(message)); %disp(p1);
Symbol_Error_Trial(k)=p1;

% Cluster Variance Calculation
ideal_symbols = [-3, -1, 1, 3];
cluster_variances = zeros(size(ideal_symbols));
for i = 1:length(ideal_symbols)
    symbol_indices = find(mprime == ideal_symbols(i));
    cluster_variances(i) = var(LWPfilteredSignal(symbol_indices));
end
%disp('Cluster Variances for each symbol:');
%disp(cluster_variances);

end
Symbol_Error_per_offset(j)=mean(Symbol_Error_Trial);
end
% % Constellation Diagram
% figure; % Create a new figure for the constellation diagram
% scatter(real(LWPfilteredSignal(1:length(LWPfilteredSignal))), zeros(1, length(LWPfilteredSignal)), 'filled');
% xlabel('Real Part');
% ylabel('Imaginary Part');
% title('Constellation Diagram ');
% grid on;
% 
% % Figures
% figure; % Create a new figure for these plots
% subplot(511)
% plot(t,s);
% title('Original Signal');
% subplot(512)
% plot(t,rsc);
% title('Carrier Added Signal');
% subplot(513)
% plot(t,noisy_rsc);
% title('Noisy RSC Signal');
% subplot(514)
% plot(t,q);
% title('Squared Signal');
% subplot(515)
% plot(t,rp);
% title('BPF Filtered Signal');
% 
% % Spectrum Plots
% figure; % Create a new figure for each spectrum plot
% plotspec(s,Ts,"Original Signal");
% figure;
% plotspec(rsc,Ts,"Suppressed Carrier");
% figure;
% plotspec(q,Ts, "Squared Signal");
% figure;
% plotspec(rp,Ts, "BPF Filtered Signal");
% 
% % More Figures
% figure; % Create a new figure for these plots
% subplot(511)
% plot(t,s);
% title('Original Signal');
% subplot(512)
% plot(t,rsc);
% title('Carrier Added Signal');
% subplot(513)
% plot(t,noisy_rsc);
% title('Noisy RSC Signal');
% subplot(514)
% plot(t,demod_signal);
% title('Demodulated Signal');
% subplot(515)
% plot(t_adjusted,LWPfilteredSignal);
% title('LWP Filtered Signal');

figure;
plot(fofset,Symbol_Error_per_offset,"Marker","o",LineWidth=1.0);
xlabel('Frequency Offset (Hz)');
ylabel('Normalized Symbol Error');
title('Normalized Symbol Error per Frequency Offset');


% Function to plot the spectrum is modified for adding labels to plots.
function plotspec(x, Ts, figureTitle)
    N = length(x);                           % length of the signal x
    t = Ts * (1:N);                          % define a time vector
    ssf = (ceil(-N/2):ceil(N/2)-1) / (Ts*N); % frequency vector
    fx = fft(x(1:N));                        % do DFT/FFT
    fxs = fftshift(fx);                      % shift it for plotting
    fxs = fxs * Ts;
    subplot(2,1,1), plot(t, x)               % plot the waveform
    xlabel('seconds'); ylabel('amplitude')   % label the axes
    title('Amplitude of ' + figureTitle);
    subplot(2,1,2), plot(ssf, abs(fxs))      % plot magnitude spectrum
    xlabel('frequency'); ylabel('magnitude') % label the axes
    title('Frequency Response of ' + figureTitle);
end