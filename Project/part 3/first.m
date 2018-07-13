clear all
% Number of iterations for each SNR
N_iter=100;

% Codeword length of LDPC
n_ldpc=64800;

% LDPC code rate
R1=2/5;
R2=9/10;
R3=2/3;
R4=8/9;

% LDPC message length
k_ldpc1=25920;
k_ldpc2=58320;
k_ldpc3=43200;
k_ldpc4=57600;

SNR= [];
SNR(1,1:14)= [-2.4:0.1:-1.1];
SNR(1,15:29)=[-1:0.1:0.4];
SNR(1,30:41)=[0.5:0.5:6];
SNR(1,42:62)=[6.1:0.02:6.5];
SNR(1,63:83)=[6.5:0.01:6.7];
SNR(1,84:90)=[7:0.5:10];
SNR(1,91:101)=[10.3:0.03:10.6];
SNR(1,102:110)=[11:0.5:15];

% LDPC Encoding and Decoding Objects
hEnc1=comm.LDPCEncoder('ParityCheckMatrix',dvbs2ldpc(R1));
hEnc2=comm.LDPCEncoder('ParityCheckMatrix',dvbs2ldpc(R2));
hEnc3=comm.LDPCEncoder('ParityCheckMatrix',dvbs2ldpc(R3));
hEnc4=comm.LDPCEncoder('ParityCheckMatrix',dvbs2ldpc(R4));
hDec1=comm.LDPCDecoder('ParityCheckMatrix',dvbs2ldpc(R1));
hDec2=comm.LDPCDecoder('ParityCheckMatrix',dvbs2ldpc(R2));
hDec3=comm.LDPCDecoder('ParityCheckMatrix',dvbs2ldpc(R3));
hDec4=comm.LDPCDecoder('ParityCheckMatrix',dvbs2ldpc(R4));

% Modulator object
M1=4;
M2=8;
hMod1=comm.PSKModulator(M1,'BitInput',true);
hMod2=comm.PSKModulator(M2,'BitInput',true);
n=1;
for gamma=SNR
    [R1, gamma, M1]
    hChan = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (SNR)','SNR',gamma);
    hDemod1=comm.PSKDemodulator(M1, 'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'Variance',1/10^(hChan.SNR/10));
    hError=comm.ErrorRate;
    EE=[];
    for counter=1:N_iter
        data=logical(randi([0,1],k_ldpc1,1));
        encodedDate=step(hEnc1,data);
        modSignal=step(hMod1,encodedDate);
        receivedSignal=step(hChan,modSignal);
        demodSignal=step(hDemod1,receivedSignal);
        receivedBits=step(hDec1,demodSignal);
        errorStats=step(hError,data,receivedBits);
        EE(1,counter)=errorStats(1);
    end
    E1(1,n)=mean(EE(1,:));
    n=n+1;
end
n=1;
for gamma=SNR
    [R2, gamma, M1]
    hChan = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (SNR)','SNR',gamma);
    hDemod1=comm.PSKDemodulator(M1, 'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'Variance',1/10^(hChan.SNR/10));
    hError=comm.ErrorRate
    EE=[];
    for counter=1:N_iter
        data=logical(randi([0,1],k_ldpc2,1));
        encodedDate=step(hEnc2,data);
        modSignal=step(hMod1,encodedDate);
        receivedSignal=step(hChan,modSignal);
        demodSignal=step(hDemod1,receivedSignal);
        receivedBits=step(hDec2,demodSignal);
        errorStats=step(hError,data,receivedBits);
        EE(1,counter)=errorStats(1);
    end
    E2(1,n)=mean(EE(1,:));
    n=n+1;
end
n=1;
for gamma=SNR
    [R3, gamma, M2]
    hChan = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (SNR)','SNR',gamma);
    hDemod2=comm.PSKDemodulator(M2, 'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'Variance',1/10^(hChan.SNR/10));
    hError=comm.ErrorRate;
    EE=[];
    for counter=1:N_iter
        data=logical(randi([0,1],k_ldpc3,1));
        encodedDate=step(hEnc3,data);
        modSignal=step(hMod2,encodedDate);
        receivedSignal=step(hChan,modSignal);
        demodSignal=step(hDemod2,receivedSignal);
        receivedBits=step(hDec3,demodSignal);
        errorStats=step(hError,data,receivedBits);
        EE(1,counter)=errorStats(1);
    end
    E3(1,n)=mean(EE(1,:));
    n=n+1;
end
n=1;
for gamma=SNR
    [R4, gamma, M2]
    hChan = comm.AWGNChannel(...
        'NoiseMethod','Signal to noise ratio (SNR)','SNR',gamma);
    hDemod2=comm.PSKDemodulator(M2, 'BitOutput',true,...
        'DecisionMethod','Approximate log-likelihood ratio',...
        'Variance',1/10^(hChan.SNR/10));
    hError=comm.ErrorRate;
    EE=[];
    for counter=1:N_iter
        data=logical(randi([0,1],k_ldpc4,1));
        encodedDate=step(hEnc4,data);
        modSignal=step(hMod2,encodedDate);
        receivedSignal=step(hChan,modSignal);
        demodSignal=step(hDemod2,receivedSignal);
        receivedBits=step(hDec4,demodSignal);
        errorStats=step(hError,data,receivedBits);
        EE(1,counter)=errorStats(1);
    end
    E4(1,n)=mean(EE(1,:));
    n=n+1;
end
semilogy(SNR,E1);
hold on
semilogy(SNR,E2);
hold on
semilogy(SNR,E3);
hold on
semilogy(SNR,E4);
grid on;
title('Example of Communication System Simulator');
ylabel('BER');
xlabel('EbNo Ratio');
legend('QPSK 2/5','QPSK 9/10','8PSK 2/3','8PSK 8/9')

