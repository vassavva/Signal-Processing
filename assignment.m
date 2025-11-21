%Task 1
[x,fs] = audioread('Vasiliki Savva.wav');
Ts = 1/fs;                      %sample period
t = 0:Ts:(length(x)-1)*Ts;      %generate discrete time values (nTs)
figure; plot(t,x);              %plot 
xlabel('t'); ylabel('x'); title('Received AM Signal');
grid on

