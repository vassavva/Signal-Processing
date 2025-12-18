%Task 1
clc;
close all;
clear all;

[x,fs] = audioread('Vasiliki Savva.wav');   %extract signal from audio
Ts = 1/fs;                                  %sampling period
N = length(x);                              %number of samples
t = 0:Ts:(N-1)*Ts;                          %generate discrete time values (nTs)
figure; plot(t,x);                          %plot time-domain
xlabel('time(sec)'); ylabel('Amplitude'); title('Time Domain');