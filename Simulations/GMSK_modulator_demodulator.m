%**************************************************************************
%
%                   written by : Abhishek Bhatta
%
% This code shows how to transmit and receive a GMSK modulated signal in
% the presence of a particular channel with additive white gaussian noise
% (AWGN). The comparison is shown between the SNR and BER of the equalized
% and non equalized signal. The main purpose of this implementation is to
% see how the GSM signal works in the presence of channel.
%
% To run this code in Octave download and install the communications
% package from the Octave sourceforge website and load the package with the
% following command (without the quotes):
%
%                       'pkg load communications'
%
%
%**************************************************************************

clear all;
close all;
clc;

% Input Parameters

SNR = (-10:2:20)';                            % SNR Values in dB
N_iterations = input('Number of Iterations : ');         % maximum number of iterations
eq_selection = input('Select the equalizer (1=LS,2=ZFE) : ');
no_of_taps = input('Select the number of taps in the channel (5,8,11,14 taps) : ');
BT = 0.3;                                   % BT product of the filter
Tb = 1e-6;                                  % bit duration
B = BT/Tb;                                  % Bandwidth
sps = 36;                                   % samples per symbol
Ts = Tb/sps;                                % sample period
M = 0.5;                                    % Modulation Index
fc = 900e6;                                 % Career frequency
l = 1e-10;                                  % changes based on the valur of fc

for snr_calc = 1:length(SNR)
    BER_eq = [];                            %BER initialisation
    BER_normal = [];
    for er = 1:N_iterations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Data Generation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        data = [0 0 1 0 0 1 0 1 1 1 0 0 0 0 1 0 0 0 1 0 0 1 0 1 1 1];   %Training Sequence(1 of 8 available in GSM)

        for k = 1:length(data)
            if data(k) == 0
                nrz_data(k) = -1;
            elseif data(k) == 1
                nrz_data(k) = 1;
            end 
        end
        
        gauss = gaussian_pulse_shaping_filter(BT,sps,Tb);       %Gaussian filter response
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Modulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        nrz = upsample(nrz_data, sps);              %Upsampling the input symbols
        nrz_gauss = conv(gauss, nrz);               %Filtering with Gaussian response
        nrz_int = cumsum(nrz_gauss);                %Integrating the signal
        nrz_gmsk = exp(1i*nrz_int);           %Generating the signal
        I_data = imag(nrz_gmsk);
        Q_data = real(nrz_gmsk);
        
        k = length(nrz_gauss)*l;
        t1 = l:l:k;
        tx_signal = cos(2*pi*fc*t1).*Q_data - sin(2*pi*fc*t1).*I_data;      %The signal to be transmitted

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Channel Impairements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ('h' is the Channel Impulse Response coefficients)
%   When choosing a particular value of 'h' make sure the 'eq_rx_signal'
%   command is also changed which appears after cancelling the channel effects.
    if no_of_taps == 5
          h = [0.0007,0.2605,0.9268,0.2605,0.0007]';              
%           eq_rx_signal = rx_signal(3:end-2);
%           eq_rx_signal = rx_signal(5:end-4);          When Equalizing
 
    elseif no_of_taps == 8

          h = [0.0007,0.0315,0.2605,0.7057,0.9268,0.7057,0.2605,0.0315]';                 
%           eq_rx_signal = rx_signal(5:end-3);
%           eq_rx_signal = rx_signal(9:end-6);          When Equalizing
  
    elseif no_of_taps == 11

          h = [0.0007,0.0108,0.0756,0.2605,0.5542,0.8272,0.9268,0.8272,0.5542,0.2605,0.0756]';        
%           eq_rx_signal = rx_signal(7:end-4);
%           eq_rx_signal = rx_signal(13:end-8);         When Equalizing


    elseif no_of_taps == 14
          h = [0.0007,0.0061,0.0315,0.1076,0.2605,0.4789,0.7057,0.8692,0.9268,0.8692,0.7057,0.4789,0.2605,0.1076]';       
%           eq_rx_signal = rx_signal(9:end-5);  
%           eq_rx_signal = rx_signal(16:end-11);        When Equalizing 

    else
        fprintf('Error: Select the number of taps from the mentioned values')
    end
        
        rx_signal = conv(h,tx_signal);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% AWGN %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        awgn_noise = awgn(rx_signal, SNR(snr_calc));           % adds white gaussian noise to the input signal
        rx_signal = rx_signal + awgn_noise;
        
        
        LP = fir1(20,0.05);                               % Low pass filter to be used for filtering the I and the Q signal.
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Demodulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %   Without Equalizing the channel

        k_normal = length(rx_signal)*l;
        t_normal = l:l:k_normal;
        I_rx_normal = rx_signal.*(-sin(2*pi*fc*t_normal));
        Q_rx_normal = rx_signal.*(cos(2*pi*fc*t_normal));

        I_rx_normal = conv(LP,I_rx_normal);
        Q_rx_normal = conv(LP,Q_rx_normal);

        phase_normal = atan2(I_rx_normal,Q_rx_normal);
        deri_normal = diff(phase_normal);
        downsampled_deri_normal = downsample(deri_normal,sps);

        %   Converting to bit sequence
        for k = 1:length(downsampled_deri_normal),
            if downsampled_deri_normal(k)>0
                data_out_normal(k) = 1;
            elseif downsampled_deri_normal(k)<0
                data_out_normal(k) = 0;
            end
    
        end


        %   Bit Error Rate Calculations
        [n_normal,r_normal]=biterr(data,data_out_normal(3:length(data_out_normal)-3));
        BER_normal = [BER_normal r_normal];
        
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Channel Estimation and Equalization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if eq_selection == 1
    % Least Square Estimation
        h2 = length(h);                             % Ensuring the length of the Toeplitz matrix
        X = toepmat(tx_signal,h2);                  %creating the toeplitz matrix of the received data
        X1 = conj(X);
        Xt = X1';                                   % Hetmitian transpose matrux is the transpose of the complex conjugate of a matrix
        inverse = inv(Xt*X);
        multi = (Xt'*inverse);
        w = (rx_signal*multi)'./2;                     % Channel Estimation values for each iteration
    elseif eq_selection == 2
    % Zero Forcing Equalizer
        N1 = 300;
        Ord = length(h);
        Lambda = 0.98;
        delta = 0.001;
        P = delta*eye(Ord);
        w = zeros(Ord,1);
        for n = Ord:N1
            u = tx_signal(n:-1:n-Ord+1);
            pin = P*u';
            k = Lambda + u*pin;
            K = pin/k;
            e(n) = rx_signal(n) - w'*u';
            w = w + K *e(n);
            PPrime = K*pin';
            P = (P-PPrime)/Lambda;
            w_err(n) = norm(h-w);
        end
    end
    
        

       rx_signal = conv(1/w,rx_signal);         % cancelling the channel efects
       
%        eq_rx_signal = rx_signal(5:end-4);       % Put the value according to the selected CIR 'h'
    if no_of_taps == 5
%           eq_rx_signal = rx_signal(3:end-2);
          eq_rx_signal = rx_signal(5:end-4);          %When Equalizing
 
    elseif no_of_taps == 8

%           h = [0.0007,0.0315,0.2605,0.7057,0.9268,0.7057,0.2605,0.0315]';                 
%           eq_rx_signal = rx_signal(5:end-3);
          eq_rx_signal = rx_signal(9:end-6);          %When Equalizing
  
    elseif no_of_taps == 11

%           h = [0.0007,0.0108,0.0756,0.2605,0.5542,0.8272,0.9268,0.8272,0.5542,0.2605,0.0756]';        
%           eq_rx_signal = rx_signal(7:end-4);
          eq_rx_signal = rx_signal(13:end-8);         %When Equalizing


    elseif no_of_taps == 14
%           h = [0.0007,0.0061,0.0315,0.1076,0.2605,0.4789,0.7057,0.8692,0.9268,0.8692,0.7057,0.4789,0.2605,0.1076]';       
%           eq_rx_signal = rx_signal(9:end-5);  
          eq_rx_signal = rx_signal(16:end-11);        %When Equalizing 
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Demodulation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       %    After cancelling channel effects

        I_rx_eq = eq_rx_signal.*(-sin(2*pi*fc*t1));
        Q_rx_eq = eq_rx_signal.*(cos(2*pi*fc*t1));

        I_rx_eq = conv(LP,I_rx_eq);
        Q_rx_eq = conv(LP,Q_rx_eq);

        phase_eq = atan2(I_rx_eq,Q_rx_eq);
        deri_eq = diff(phase_eq);
        downsampled_deri_eq = downsample(deri_eq,sps);

       %    Converting to bit sequence
        for k = 1:length(downsampled_deri_eq),
            if downsampled_deri_eq(k)>0
                data_out_eq(k) = 1;
            elseif downsampled_deri_eq(k)<0
                data_out_eq(k) = 0;
            end
    
        end
        
       %    Bit Error Rate Calculations
        [n_eq,r_eq]=biterr(data,data_out_eq(3:length(data_out_eq)-3));
        BER_eq = [BER_eq r_eq];
    

    end
    
    ber_total_eq(snr_calc) = sum(BER_eq)/N_iterations;
    ber_total_normal(snr_calc) = sum(BER_normal)/N_iterations;
    
end

ber_eq = [SNR ber_total_eq']                % BER vs SNR for equalized signal
ber_normal = [SNR ber_total_normal']        % BER vs SNR without equalization



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Plotting the Figures %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(1)
% semilogy(SNR,ber_total_eq,'linewidth',2);hold on
% semilogy(SNR,ber_total_normal, '-*','linewidth',2);
% set(gca,'FontSize',16)
% xlabel('SNR [dB]')
% ylabel('BER')
% legend('Equalized Output','Without Equalization')
% grid on
% hold off
% 
% figure(2);
% subplot(2,2,1);stem(data);title('Input data');xlabel('Time');ylabel('Amplitude');
% subplot(2,2,2);plot(I_data);title('Transmitted IQ data');xlabel('Time');ylabel('Amplitude');hold on
% plot(Q_data, 'r');xlabel('Time');ylabel('Amplitude');
% subplot(2,2,3);plot(I_rx_eq);title('Received IQ data after filtering');xlabel('Time');ylabel('Amplitude');hold on
% plot(Q_rx_eq, 'r');xlabel('Time');ylabel('Amplitude');
% subplot(2,2,4);stem(data_out_eq);title('Demodulated data');xlabel('Time');ylabel('Amplitude');
%  
% figure(3)
% stem(nrz_out_eq(3:length(nrz_out_eq)-3),'*');hold on
% stem(nrz_data)
%   
figure(4)
plot(h);hold on
plot(w,'r');
legend('Channel Impairements','Estimated Channel')
% 
% figure(5)
% scatter(I_data,Q_data);hold on
% scatter(I_rx_eq,Q_rx_eq)
% legend('Tx I Q data','Rx I Q data after filtering')
% xlabel('In Phase');ylabel('Quadrature')
% 
% figure
% stem(rx_signal)
% hold on
% stem(tx_signal,'*')
% legend('Equalized Received Signal','Transmitted Signal')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
