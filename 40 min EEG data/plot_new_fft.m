[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
[da,ca]= butter(3,[38/250 40/250], 'stop');
load('CoStim_Result_Subjectfdz2_Sess2_Stim2000mA.mat')

START_ses2_CH4 = 18910;
END_ses2_CH4 =  START_ses2_CH4 + 1199999;

figure(987);
for i = 1:1:8
    DATA_SES_1=filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(i,:)));
    subplot(4,2,i)
    plot(DATA_SES_1/10000000)
    xlabel('samples')
    ylabel('EEG/ mV')
    xlim([1 1250000])
end


DATA_CH4=filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(4,START_ses2_CH4:END_ses2_CH4)));

% FFT settings

f_samp_b = 500;
k = 1;
for i=START_ses2_CH4:150000:1068910
low = i;
high = i+149999;
    
DATA_5min = filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(4,low:high)));

eeg_b =DATA_5min/10000000;

epoch_duration = 1; n_fft = 2^15; windows = f_samp_b * epoch_duration; overlap = 0.1 * f_samp_b;

eeg_b_to_analyse = eeg_b(1, 40*f_samp_b:end);

[pxx_b, f_b] = pwelch(eeg_b_to_analyse,windows,overlap,n_fft,f_samp_b);

figure(3);

subplot(4,2,k)
semilogx(f_b,10*log10(pxx_b),'LineWidth',2,'Color','b');
xlabel('f (Hz)')
ylabel('|P1(f)| (dB)')
xlim([1 100])

k = k+1;
end 



