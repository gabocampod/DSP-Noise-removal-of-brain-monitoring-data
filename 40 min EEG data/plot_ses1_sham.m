[b,a]=butter(3,0.5/250,'high');
[d,c]=butter(3,50/250,'low');
[da,ca]= butter(3,[38/250 40/250], 'stop');
load('CoStim_Result_Subjectfdz1_Sess1_Stim2000mA.mat')


sham_tacs_t1 = 320000;
sham_tacs_t2 = sham_tacs_t1 + 149999;

figure(1);
for i = 1:1:8
    Sham_SMA=filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(i, :)));
    subplot(4,2,i)
    plot(Sham_SMA)
end

% FFT settings

f_samp_b = 500;
k = 1;

for i=1:150000:1050001
low = i;
high = i+149999;
    
DATA_5min_sHAM = filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(4,low:high)));

eeg_b =DATA_5min_sHAM;

epoch_duration = 1; n_fft = 2^15; windows = f_samp_b * epoch_duration; overlap = 0.1 * f_samp_b;

eeg_b_to_analyse = eeg_b(1, 40*f_samp_b:end);

[pxx_b, f_b] = pwelch(eeg_b_to_analyse,windows,overlap,n_fft,f_samp_b);

figure(2);
subplot(4,2,k)
semilogx(f_b,10*log10(pxx_b),'LineWidth',2,'Color','b'); hold all

k = k+1;
end 





