clear;
% Convolutional Encoder with Puncturing at the code rate of 1/2 
hConvEnc1 = comm.ConvolutionalEncoder(poly2trellis(7,[133 171]));
hConvEnc1.PuncturePatternSource = 'Property';
hConvEnc1.PuncturePattern = [1;1];

% Convolutional Encoder with Puncturing at the code rate of 2/3 
hConvEnc2 = comm.ConvolutionalEncoder(poly2trellis(7,[133 171]));
hConvEnc2.PuncturePatternSource = 'Property';
hConvEnc2.PuncturePattern = [1;1;0;1];

% Convolutional Encoder with Puncturing at the code rate of 5/6
hConvEnc3 = comm.ConvolutionalEncoder(poly2trellis(7,[133 171]));
hConvEnc3.PuncturePatternSource = 'Property';
hConvEnc3.PuncturePattern = [1;1;0;1;1;0;0;1;1;0];

% Viterbi Decoder and the punctured pattern at the code rate of 1/2
hVitDec1 = comm.ViterbiDecoder(poly2trellis(7,[133 171]),'InputFormat','Hard');
hVitDec1.PuncturePatternSource = 'Property';
hVitDec1.PuncturePattern = hConvEnc1.PuncturePattern; 
hVitDec1.TracebackDepth = 96;
hErrorCalc1 = comm.ErrorRate('ReceiveDelay',hVitDec1.TracebackDepth);

% Viterbi Decoder and the punctured pattern at the code rate of 2/3
hVitDec2 = comm.ViterbiDecoder(poly2trellis(7,[133 171]),'InputFormat','Hard');
hVitDec2.PuncturePatternSource = 'Property';
hVitDec2.PuncturePattern = hConvEnc2.PuncturePattern;
hVitDec2.TracebackDepth = 96;
hErrorCalc2 = comm.ErrorRate('ReceiveDelay',hVitDec2.TracebackDepth);

% Viterbi Decoder and the punctured pattern at the code rate of 5/6
hVitDec3 = comm.ViterbiDecoder(poly2trellis(7,[133 171]),'InputFormat','Hard');
hVitDec3.PuncturePatternSource = 'Property';
hVitDec3.PuncturePattern = hConvEnc3.PuncturePattern;
hVitDec3.TracebackDepth = 96;
hErrorCalc3 = comm.ErrorRate('ReceiveDelay',hVitDec3.TracebackDepth);

% Modulator and Demodulator
hMod = comm.QPSKModulator('BitInput',true);
hDemod = comm.QPSKDemodulator('BitOutput',true);

% AWGN channels
hChan1 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)',...
    'SignalPower',1,'SamplesPerSymbol',1,'BitsPerSymbol',2);
hChan2 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)',...
    'SignalPower',1,'SamplesPerSymbol',1,'BitsPerSymbol',2);
hChan3 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)',...
    'SignalPower',1,'SamplesPerSymbol',1,'BitsPerSymbol',2);


% Define the SNR range needed to ensure BER range at least from 10-5 to 10-2 can be obtained 
EbNoEncoderInput = 0:.5:15;

% Eb/No ration incorporating puncturing rates
EbNoEncoderOutput1 = EbNoEncoderInput + 10*log10(1/2);
EbNoEncoderOutput2 = EbNoEncoderInput + 10*log10(2/3);
EbNoEncoderOutput3 = EbNoEncoderInput + 10*log10(5/6);
frameLength = 12000;
targetErrors = 100;
maxNumTransmissions = 5e6;

% Allocate memory to store results to keep track on the simulation status
BERVec1 = zeros(3,length(EbNoEncoderInput));
BERVec2 = zeros(3,length(EbNoEncoderInput));
BERVec3 = zeros(3,length(EbNoEncoderInput));


%% Simulation Loop
for n = 1:length(EbNoEncoderInput)
    
    % reset all Matlab handler to default
    reset(hErrorCalc1)
    reset(hConvEnc1)
    reset(hVitDec1)
    reset(hErrorCalc2)
    reset(hConvEnc2)
    reset(hVitDec2)  
    reset(hErrorCalc3)
    reset(hConvEnc3)
    reset(hVitDec3)
   
    % generate AWGN Channel for each Eb/No value
    hChan1.EbNo = EbNoEncoderOutput1(n);
    hChan2.EbNo = EbNoEncoderOutput2(n);
    hChan3.EbNo = EbNoEncoderOutput3(n);
 
    
    while (((BERVec1(2,n) < targetErrors)...
        ||(BERVec2(2,n) < targetErrors)||(BERVec3(2,n) < targetErrors))...
        && (BERVec3(3,n) < maxNumTransmissions))
        % Generate binary frames of size specified by the frameLength
        % variable
        data = randi([0 1],frameLength,1);
        
        % Convolutionally encode the interleaved data
        encData1 = step(hConvEnc1,data);
        encData2 = step(hConvEnc2,data);
        encData3 = step(hConvEnc3,data);
        
        % Modulate the encoded data into a QPSK symbol
        modData1 = step(hMod,encData1);
        modData2 = step(hMod,encData2);
        modData3 = step(hMod,encData3);

        
        % Add a Rayleigh channel
        Rayleigh1 = sqrt(0.5*((randn(1,length(modData1))).^2 + (randn(1,length(modData1))).^2)); 
        Rayleigh2 = sqrt(0.5*((randn(1,length(modData2))).^2 + (randn(1,length(modData2))).^2));
        Rayleigh3 = sqrt(0.5*((randn(1,length(modData3))).^2 + (randn(1,length(modData3))).^2)); 
        
        ChannalOut1 = Rayleigh1'.*modData1;
        ChannalOut2 = Rayleigh2'.*modData2;
        ChannalOut3 = Rayleigh3'.*modData3;
        
        % Pass the modulated signal through an AWGN channel
        channelOutput1 = step(hChan1,ChannalOut1);
        channelOutput2 = step(hChan2,ChannalOut2);
        channelOutput3 = step(hChan3,ChannalOut3);

        
        % Demodulation QPSK symbols back to encoded data
        DemodData1 = step(hDemod,channelOutput1);
        DemodData2 = step(hDemod,channelOutput2);
        DemodData3 = step(hDemod,channelOutput3);

        
        % Input to the hard-decision Viterbi decoder
        decData1 = step(hVitDec1,DemodData1);
        decData2 = step(hVitDec2,DemodData2);
        decData3 = step(hVitDec3,DemodData3);
        
        % Compute the errors for each data with the original data
        BERVec1(:,n) = step(hErrorCalc1,data,decData1);
        BERVec2(:,n) = step(hErrorCalc2,data,decData2);
        BERVec3(:,n) = step(hErrorCalc3,data,decData3);

    end
end
% save the data
BER_QPSK1 = BERVec1(1,:);
BER_QPSK2 = BERVec2(1,:);
BER_QPSK3 = BERVec3(1,:);
save('QPSK_sim3.mat','EbNoEncoderInput','BER_QPSK1','BER_QPSK2','BER_QPSK3');

%% Plot Results
load('QPSK_sim3.mat');
hold on;
semilogy(EbNoEncoderInput,BER_QPSK1);
hold on
semilogy(EbNoEncoderInput,BER_QPSK2);
hold on
semilogy(EbNoEncoderInput,BER_QPSK3);
grid on;
title('Example of Communication System Simulator');
ylabel('BER');
xlabel('EbNo Ratio');
legend('1/2','2/3','5/6')