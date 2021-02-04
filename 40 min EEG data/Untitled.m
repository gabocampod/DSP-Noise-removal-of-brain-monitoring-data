figure(122);
semilogx(f_b_sham,10*log10(pxx_b_sham)-60,'LineWidth',2,'Color','b'); hold all
semilogx(f_b,10*log10(pxx_b)-60,'LineWidth',2,'Color','r'); hold all
xlabel('f(Hz)')
ylabel('Magnitude (dB)')
xlim([1 100])