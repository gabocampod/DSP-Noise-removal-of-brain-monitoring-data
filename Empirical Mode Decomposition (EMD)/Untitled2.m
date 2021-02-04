figure(10);
for m=1:50
  subplot(10,5,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  %ylim([-1 1]);
end;