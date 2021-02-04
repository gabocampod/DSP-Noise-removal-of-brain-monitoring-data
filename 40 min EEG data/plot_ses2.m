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
    plot(DATA_SES_1)
end
title('DATA SES 2')

DATA_CH4=filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(4,START_ses2_CH4:END_ses2_CH4)));

figure(398);
plot(DATA_CH4)

Fs = 500;            %FFT of fILTERED DATA                          
L = 1200000;            
Y = fft(DATA_CH4);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L;
figure(2);
plot(f,P1, 'k') 
xlabel('f(Hz)')
ylabel('|P1(f)|')



