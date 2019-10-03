clc
clear 
close all    

%% Parameters

sec = 60;  % data plot in seconds
    
%     load('C:\Users\Asus\Downloads\Biomed Project\SP Cup 2014 Dataset\competition_data\Training_data\DATA_01_TYPE01.mat')

for i=2:9 
    load(['C:\Users\Asus\Downloads\Biomed Project\SP Cup 2014 Dataset\competition_data\Training_data\DATA_0',num2str(i),'_TYPE0',num2str(2),'.mat'])
    
    ecg = sig(1,:);   %ecg signal
    ppg1= sig(2,:);   %ppg signal
    ppg2= sig(3,:);   %ppg signal
    ax = sig(4,:);  %acceleration signal along x axis
    ay = sig(5,:);  %acceleration signal along y axis
    az = sig(6,:);  %acceleration signal along z axis

    fs = 125;

    subplot(211),plot(ppg1(1:sec*fs))
    title(['Person-',num2str(i),'  ppg signal channel-1'])
    subplot(212),plot(ppg2(1:sec*fs))
    title(['Person-',num2str(i),'  ppg signal channel-2'])
    pause(2)
    
end

