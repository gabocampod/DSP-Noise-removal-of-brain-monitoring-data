[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');

%top lim = 1218870
load('CoStim_Result_Subjectfdz2_Sess2_Stim2000mA.mat')
tacs_t1 = 168920;        
tacs_t2 = tacs_t1 + 149999; %5 MINUTES (5*60*500)

                                                             
[ReconSignal_SUBJECT] = SMA_ERP2(RAW_SUBJECT_DATA(1,tacs_t1: tacs_t2),40);  %%CHANGE FREQ
FILTERED_FREE_EEG = filtfilt(b,a,filtfilt(d,c,ReconSignal_SUBJECT(1,:)));

figure(1, plot)