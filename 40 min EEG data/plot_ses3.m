[d,c]=butter(3,50/250,'low');
[b,a]=butter(3,0.5/250,'high');
[da,ca]= butter(3,[38/250 40/250], 'stop');
load('CoStim_Result_Subjectfdz3_Sess3_Stim2000mA.mat')

start_ses3_ch4 = 16019;
end_ses3_ch4 = start_ses3_ch4 + 1199999;

figure(1);
for i = 1:1:8
    DATA_SES_3=filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(i,:)));
    subplot(4,2,i)
    plot(DATA_SES_3)
end
title('DATA SES 3')

DATA_CH4 = filtfilt(b,a,filtfilt(d,c,CoStim_Result.EEG.Data(4,start_ses3_ch4:end_ses3_ch4)));
%DATA_CH4 = filtfilt(da,ca,DATA_CH4);

figure(3);
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
title('FFT CHAN 4')






