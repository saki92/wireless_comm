classdef SumRate
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Nt = 2; %no of transmit antennas
        K = 2; %no of receivers
        Q = 1; %no of antennas per receiver
        M; %no of error samples
        MV; %no of error samples for rate validation
        sigma = 1; %noise covariance
        SNR = 5:2:35; %dB
        PtL = 10.^(SNR/10); %total transmit power
        noofit; %no of validation iteration
        error_var = 0.3; %alpha factor
    end
    
    methods
    end
    
end

