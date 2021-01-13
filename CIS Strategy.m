addpath('D:\Carleton\Courses\Signal Processing\articles\codes');
[Data,Srate] = audioread('df1_n2H.wav');
%%%%%%%%%%% checking if the signal is mono or stereo %%%%%%%%%%% 
[~, c]=size(Data);
if c>1
    MonoSig = sum(Data,2);
else
    MonoSig = Data;
end

figure(1), clf
plot(MonoSig);
title('Input Waveform');
set(gca,'ylim',[-1.2 1.2])
xlabel('Sample Number');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%changing sampling rate to 16KHz %%%%%%%%%%%%%%

[P,Q] = rat(16000/Srate);
NewSig = resample(MonoSig,P,Q);
NewSig = NewSig(:,1);
fs = P/Q*Srate;

figure(2), clf
plot(NewSig);
title('Input Waveform with 16KHz Srate');
set(gca,'ylim',[-1.2 1.2])
xlabel('Sample Number');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%% Pre-Emphesizing %%%%%%%%%%%%%%%%%%%%

fc    = 1200;
w     = 2*fc/fs; % fc/fs/2 = 2*fc/fs = 1200 Hz
[b,a] = butter(2,w,'high');
SigPreEmp = filter(b,a,NewSig);
%SigPreEmp = NewSig;

figure(3), clf
plot(SigPreEmp);
title('Pre-Emphesized Input Waveform ');
set(gca,'ylim',[-1.2 1.2])
xlabel('Sample Number');
ylabel('Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%% Channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Channels = 8;
ChBand = 500;
OutputSig = 0;
for i = 1:Channels    
    % Set lower and upper corner frequency of bandpass filter
    if i == 1
        lowerband = 100;
        higherband = 500;
    else
        lowerband = 500 + (i-2)*ChBand;
        higherband = lowerband + ChBand;
    end
    
    filterOrder = 2;
    BandFiltSig = butterBandpassFilter(SigPreEmp, lowerband, higherband, fs, filterOrder);
    CenterFreq = (lowerband + higherband) / 2;
    
  %figure(10)  
   %     subplot (8,1,i)
    %    plot(BandFiltSig);
     %   suptitle('Band-Pass Filters');
      %  ylabel(['Ch',num2str(i)]);
       % set(gca,'ylim',[-0.2 0.2])
        

    
    %%%%% rectify each band %%%%%%
   RectSig = abs(BandFiltSig);
   %%%%% Envelope Extraction for each band %%%%%
   SigEnvelope = butterLowpassFilter(RectSig, 400, fs, filterOrder);
   %figure(11)  
    %    subplot (8,1,i)
     %   plot(SigEnvelope);
      %  suptitle('Envelope of each channel');
       % ylabel(['Ch',num2str(i)]);
        %set(gca,'ylim',[-0.2 0.2])
   
   % Generate cosine signal with central frequency of bandpass filters and length of rectified signal
    [r, ~] = size(RectSig);
    timeDuration = r/fs;
    time = linspace(0, timeDuration, r);
    CosSig = cos(2*pi*CenterFreq*time);
    
    %%%%% AMPLITUDE MODULATION %%%%%
    SigEnvelope = transpose(SigEnvelope);
    ModSig = SigEnvelope.*(CosSig);
%    figure(12)  
 %       subplot (8,1,i)
  %      plot(ModSig);
   %     suptitle('Modulated Signal for each channel');
    %    ylabel(['Ch',num2str(i)]);
     %   set(gca,'ylim',[-0.2 0.2])
    
    % Sum amplitude modulated signals for each channel
    OutputSig = OutputSig + ModSig;
    
    %%%%% Plotting %%%%%
    if i == 1
        
        N=length(BandFiltSig);
        sigpow= abs( fft(BandFiltSig)/N ).^2;
        hz = linspace(0,fs/2,floor(N/2)+1);

        figure()   
        plot(hz,sigpow(1:length(hz)));
        title('Output Signal of Lowest Frequency Channel');
        xlabel('Sample Number');
        ylabel('Amplitude');
        set(gca,'xlim',[0 1600])
        
        figure()
        plot(SigEnvelope);
        title('Envelope of Lowest Frequency Channel');
        xlabel('Sample Number');
        
    elseif i == Channels
        
        N=length(BandFiltSig);
        sigpow= abs( fft(BandFiltSig)/N ).^2;
        hz = linspace(0,fs/2,floor(N/2)+1);

        figure()   
        plot(hz,sigpow(1:length(hz)));
        title('highest Frequency Channel');
        xlabel('Sample Number');
        ylabel('Amplitude');
        set(gca,'xlim',[0 8000])
        
        figure()   
        plot(BandFiltSig);
        title('Output Signal of Highest Frequency Channel');
        xlabel('Sample Number');
        ylabel('Amplitude');
        
        figure()
        plot(SigEnvelope);
        title('Envelope of Highest Frequency Channel');
        xlabel('Sample Number');
    end
        
end
% Normalize the signals by the max of their absolute value
NormOutputSig = OutputSig / max(abs(OutputSig), [], 'all');
%NormInput = SigPreEmp / max(abs(SigPreEmp), [], 'all');
NormInput = NewSig / max(abs(NewSig), [], 'all');

%%%%% Plot the normalized output signal
figure()
plot(NormOutputSig);
title('Synthsized Signal');
xlabel('Sample Number');
ylabel('Amplitude');

% Plot the normalized input signal
figure()
plot(time, NormInput);
hold on;
plot(time, NormOutputSig);
title('Original vs Processed Sound');
legend('Original', 'Processed');
sound(OutputSig, fs);

%for i = 1:Channels
 %       figure(10)  
  %      subplot (12,1,i)
   %     plot(BandFiltSig);
    %    suptitle('BP filter');
     %   ylabel(['Ch',num2str(i)]);
 % end

function[y] = butterBandpassFilter(data, lowcut, highcut, fs, order)
    % Nyquist frequency
    nyq = fs/2;

    % Since the cutoff frequency cannot be equal to 1 and nyq = 8000,
    % the upper cutoff frequency must be less than 8000
    if highcut == nyq
        highcut = highcut - 0.00000000001;
    end
    
    % Normalize the frequencies by dividing by the Nyquist frequency
    lowerband = lowcut/nyq;
    higherband = highcut/nyq;
    
    % butter() returns b,a which are transfer function coefficients
    [b, a] = butter(order, [lowerband, higherband], 'bandpass');
    y = filter(b, a, data);
end
function[y] = butterLowpassFilter(data, cutoff, fs, order)
    % Nyquist frequency 
    nyq = fs/2;
    
    % Cutoff frequency cannot equal Nyquist frequency
    % so decrease slightly
    if cutoff == nyq
        cutoff = cuttoff - 0.00000000001;
    end
    
    % Normalize the cutoff frequency
    cutoffFreq = cutoff / nyq;
    
    [b, a] = butter(order, cutoffFreq, 'low');    
    y = filter(b, a, data);
end
