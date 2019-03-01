% Author: Prabhu Chandhar, Chandhar Research Labs, Chennai, India.


clc;
clear all;
n_bits = 10000;% Number of bits
% bandwidth = 20e6;
N = 64;            % FFT length
n_sc = N;          % Number of data sub-carriers
M = 100;           % Number of antennas at the base station
U = 10;            % Number of terminals
% subcarrier_spacing = bandwidth/N; 
% ofdm_symbol_period = 1/subcarrier_spacing;    % synbols symbol duration before CP
% sampling_period = 1/(N*subcarrier_spacing); 
N_taps = 1;
cp_length = 10; 
snrdB = -10:4:20;
snr = 10.^(snrdB/10);% -15:2:40;
order_m = 4; % order of modulaton 
om = log2(order_m);
n_bits = ceil(n_bits/(om*n_sc*N))*om*n_sc*N; % Number of transmitted bits  

for lp_snr=1:length(snr)
    lp_snr
for u = 1:U
    bits_de(u).array = bi2de(reshape(randint(1,n_bits),n_bits/om,om)); % Grouping of bits for mapping
    sym_m(u).array = pskmod(bits_de(u).array,2^om);                    % Symbol mapping
end

for u = 1:U
    
    k = 1;
    for i = 1:length(sym_m(u).array)/N
        sym_block(u).array(i,:) = sym_m(u).array(k:k+N-1); % Number of symbols per block
        k = k+N;
    end
end

%**************************************OFDM Symbol Generation*************************************************************

sb = size(sym_block(1).array);
for i=1:sb(1)
    for u = 1:U
        ifft_input(u).array = zeros(1,N);  
        ifft_input(u).array = sym_block(u).array(i,:); % Input symbols to the IFFT
        ifft_output_before_cp(u).array = ifft(ifft_input(u).array,N);      % OFDM Symbol before CP
        ifft_output_after_cp(u).array = [ifft_output_before_cp(u).array((N-cp_length+1:N)) ifft_output_before_cp(u).array];  % OFDM Symbol after CP
        ofdm_samples(u).array(i,:) = ifft_output_after_cp(u).array;
    end
end

% %***********************2******************* CHANNEL **********************************************************************

    for i = 1:sb(1)
        for m = 1:M
            for u = 1:U
                channel(m).antenna(u).user(i,:) = inv(sqrt(2)*N_taps)*(randn(1,N_taps) + j*randn(1,N_taps));% Channel coefficients
                h_k(m).antenna(u).user(i,:) = fft(channel(m).antenna(u).user(i,:),N); % FFT of the channel
                                  for k = 1:N 
                                      HH(k).sc(i).block(m,u) = h_k(m).antenna(u).user(i,k); % channel coefficients
                                  end

                ofdm_samples_after_channel(m).antenna(u).user(i,:) = filter(channel(m).antenna(u).user(i,:), 1, ofdm_samples(u).array(i,:)); %filtered OFDM samples after channel
              
            end
        end
    end
    
    
    %******************************************* RECEIVER **********************
    
    for i=1:sb(1)
        for m = 1:M
            ofdm_receiver_input(m).antenna(i,1:length(ofdm_samples_after_channel(m).antenna(1).user(i,:))) = 0; %initialization
            for u = 1:U
                ofdm_receiver_input(m).antenna(i,:) = ofdm_receiver_input(m).antenna(i,:) + ofdm_samples_after_channel(m).antenna(u).user(i,:);
            end
  noise = inv(sqrt(N))*inv(sqrt(snr(lp_snr)))*inv(sqrt(2))*(randn(1,length(ofdm_receiver_input(m).antenna(i,:)))+j*randn(1,length(ofdm_receiver_input(m).antenna(i,:))));
            ofdm_receiver_input_after_adding_noise(m).antenna(i,:) = ofdm_receiver_input(m).antenna(i,:) + noise;%awgn(ofdm_receiver_input(m).antenna(i,:),snr(lp_snr),'measured'); % recieved noisy signal

            fft_input(m).antenna(i).block(1:N) = ofdm_receiver_input_after_adding_noise(m).antenna(i,cp_length+1:N+cp_length);  
            fft_output(m).antenna(i).block(1:N) = fft(fft_input(m).antenna(i).block(1:N),N);% FFT of the received signal
        Y_k(i).block(m,:) = fft_output(m).antenna(i).block(1:N);
        end

            for k = 1:N
            fft_output_hat_mrc(i).block(k,1:U) = HH(k).sc(i).block'*Y_k(i).block(:,k); % estimated signal after MRC
            fft_output_hat_zf(i).block(k,1:U) = (HH(k).sc(i).block'*HH(k).sc(i).block)^(-1)*HH(k).sc(i).block'*Y_k(i).block(:,k);% estimated signal after ZF
            end            
    end

    
    for u = 1:U
        transmitted_symbol_array = [];
        received_symbol_array_mrc = [];
        received_symbol_array_zf = [];
        for i=1:sb(1)
            transmitted_symbol_array = [transmitted_symbol_array sym_block(u).array(i,:).'];
            received_symbol_array_mrc = [received_symbol_array_mrc fft_output_hat_mrc(i).block(:,u).'];
            received_symbol_array_zf = [received_symbol_array_zf fft_output_hat_zf(i).block(:,u).'];
        end
       %*********************************Demodulation **********************   
        bits_de(u).array = pskdemod(transmitted_symbol_array,2^om);
        bits_bi(u).array = de2bi(bits_de(u).array);
        bits_i(u).array = reshape(bits_bi(u).array,1,n_bits);
        
        bits_dde_mrc(u).array = pskdemod(received_symbol_array_mrc,2^om);
        bits_dde_zf(u).array = pskdemod(received_symbol_array_zf,2^om);
        
        
        bits_dbi_mrc(u).array = de2bi(bits_dde_mrc(u).array);
        bits_o_mrc(u).array = reshape(bits_dbi_mrc(u).array,1,n_bits);
        
        bits_dbi_zf(u).array = de2bi(bits_dde_zf(u).array);
        bits_o_zf(u).array = reshape(bits_dbi_zf(u).array,1,n_bits);
        %************************************* SER Calculation ******************************************************************
        % err=bits_de-bits_dde;
        err_mrc(u).array = bits_i(u).array-bits_o_mrc(u).array;
        err_zf(u).array = bits_i(u).array-bits_o_zf(u).array;
        
        
        e_mrc(u).array = 0;
        for lp = 1:length(err_mrc(u).array)
            if err_mrc(u).array(lp)~=0
                e_mrc(u).array = e_mrc(u).array+1;
            end
            
        end
        
        BER_mrc(u).array(lp_snr) = [e_mrc(u).array/length(err_mrc(u).array)];
        
        e_zf(u).array = 0;
        for lp = 1:length(err_zf(u).array)
            if err_zf(u).array(lp)~=0
                e_zf(u).array = e_zf(u).array+1;
            end
            
        end
        
        BER_zf(u).array(lp_snr) = [e_zf(u).array/length(err_zf(u).array)];
    end
end

figure(20);semilogy(snrdB,BER_mrc(1).array,'b');grid ON; hold on;
figure(20);semilogy(snrdB,BER_zf(1).array,'k');grid ON;
xlabel('SNR');ylabel('BER');
