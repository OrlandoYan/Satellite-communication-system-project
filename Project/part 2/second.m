clear;
% Reed Solomon Encoder and Decoder
Polynomial = rsgenpoly(255,239,[],0);
hRSEnc = comm.RSEncoder(255,239,Polynomial,188);
hRSDec = comm.RSDecoder(255,239,Polynomial,188);

% Convolutional interleaver and deinterleaver
nrows = 12; slope = 17;
D = nrows*(nrows - 1)*slope;
hInt = comm.ConvolutionalInterleaver('NumRegisters',nrows,'RegisterLengthStep',slope);
hDeint = comm.ConvolutionalDeinterleaver('NumRegisters',nrows,'RegisterLengthStep',slope);

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
hVitDec1.TerminationMethod =  'Truncated';
hErrorCalc1 = comm.ErrorRate;

% Viterbi Decoder and the punctured pattern at the code rate of 2/3
hVitDec2 = comm.ViterbiDecoder(poly2trellis(7,[133 171]),'InputFormat','Hard');
hVitDec2.PuncturePatternSource = 'Property';
hVitDec2.PuncturePattern = hConvEnc2.PuncturePattern;
hVitDec2.TerminationMethod =  'Truncated';
hErrorCalc2 = comm.ErrorRate;

% Viterbi Decoder and the punctured pattern at the code rate of 5/6
hVitDec3 = comm.ViterbiDecoder(poly2trellis(7,[133 171]),'InputFormat','Hard');
hVitDec3.PuncturePatternSource = 'Property';
hVitDec3.PuncturePattern = hConvEnc3.PuncturePattern;
hVitDec3.TerminationMethod =  'Truncated';
hErrorCalc3 = comm.ErrorRate;

% Calculate the errorrate without the RS outer code 
hErrorCalc4 = comm.ErrorRate;

% Modulator and Demodulator
hMod = comm.QPSKModulator('BitInput',true);
hDemod = comm.QPSKDemodulator('BitOutput',true);

% AWGN channels
hChan1 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)',...
    'BitsPerSymbol',2);
hChan2 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)',...
    'BitsPerSymbol',2);
hChan3 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)',...
    'BitsPerSymbol',2);
hChan4 = comm.AWGNChannel('NoiseMethod','Signal to noise ratio (Eb/No)',...
    'BitsPerSymbol',2);

% Define the SNR range needed to ensure BER range at least from 10-5 to 10-2 can be obtained 
EbNoEncoderInput = -4:.5:22;

% Eb/No ration incorporating puncturing rates
EbNoEncoderOutput1 = EbNoEncoderInput + 10*log10(1/2) ;
EbNoEncoderOutput2 = EbNoEncoderInput + 10*log10(2/3) ;
EbNoEncoderOutput3 = EbNoEncoderInput + 10*log10(5/6) ;
EbNoEncoderOutput4 = EbNoEncoderInput + 10*log10(204/188);

frameLength = 188;
targetErrors = 100;
maxNumTransmissions = 5e6;

% Allocate memory to store results to keep track on the simulation status
BERVec1 = zeros(3,length(EbNoEncoderInput)); 
BERVec2 = zeros(3,length(EbNoEncoderInput)); 
BERVec3 = zeros(3,length(EbNoEncoderInput)); 
BERVec4 = zeros(3,length(EbNoEncoderInput)); 
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
    reset(hErrorCalc4)
    reset(hRSEnc)
    reset(hRSDec)
    release(hInt)
    release(hDeint)
    
   % generate AWGN Channel for each Eb/No value
    hChan1.EbNo = EbNoEncoderOutput1(n);
    hChan2.EbNo = EbNoEncoderOutput2(n);
    hChan3.EbNo = EbNoEncoderOutput3(n);
    hChan4.EbNo = EbNoEncoderOutput3(n);
    while (((BERVec1(2,n) < targetErrors)||(BERVec2(2,n) < targetErrors)||...
        (BERVec3(2,n) < targetErrors)||(BERVec4(2,n) < targetErrors)) &&... 
        (BERVec1(3,n) < maxNumTransmissions))
        % Generate binary frames of size specified by the frameLength
        % variable
        data = randi([0 1],frameLength,8);
        
        % Reed Solomon encode
        RSencoder = step(hRSEnc,bi2de(data));

        % Interleaveing process
        dataInt1 = [RSencoder;zeros(2244,1)];
        dataInt2 = step(hInt,dataInt1);
        
        % Convert output to binary bits
        dataIntbit1 = de2bi(dataInt2);
        dataIntbit2 = dataIntbit1';
        dataIntbit3 = dataIntbit2(:);
        
        % Convolutionally encode the interleaved data
        encData1 = step(hConvEnc1,dataIntbit3);
        encData2 = step(hConvEnc2,dataIntbit3);
        convert = [dataIntbit3;zeros(6,1)];
        encData3 = step(hConvEnc3,convert);
       
        % Modulate the encoded data into a QPSK symbol
        modData1 = step(hMod,encData1);
        modData2 = step(hMod,encData2);
        modData3 = step(hMod,encData3);
        
        % Modulate the uncoded data into a QPSK symbol
        modData4 = step(hMod,data(:));
        
        % Pass the modulated signal through an AWGN channel
        channelOutput1 = step(hChan1,modData1);
        channelOutput2 = step(hChan2,modData2);
        channelOutput3 = step(hChan3,modData3);
        channelOutput4 = step(hChan4,modData4);
        
        % Demodulation QPSK symbols back to encoded data
        DemodData1 = step(hDemod,channelOutput1);
        DemodData2 = step(hDemod,channelOutput2);
        DemodData3 = step(hDemod,channelOutput3);
        DemodData4 = step(hDemod,channelOutput4);
        
        % Input to the hard-decision Viterbi decoder
        decData1 = step(hVitDec1,DemodData1);
        decData2 = step(hVitDec2,DemodData2);
        decData3_1 = step(hVitDec3,DemodData3);
        decData3 = decData3_1(1:19584);
       
        % Convert back to original form
        decDatao1 = reshape(decData1,[8,2448])';
        decDatao2 = reshape(decData2,[8,2448])';
        decDatao3 = reshape(decData3,[8,2448])';
      
        % Deinterleaver process
        decDataInt1 = step(hDeint,bi2de(decDatao1));
        decDataInt2 = step(hDeint,bi2de(decDatao2));
        decDataInt3 = step(hDeint,bi2de(decDatao3));
        
        decDataInto1 = decDataInt1(2245:2448);
        decDataInto2 = decDataInt2(2245:2448);
        decDataInto3 = decDataInt3(2245:2448);
        
        % Reed Solomon decode
        deDataRS1 = step(hRSDec,decDataInto1);
        deDataRS2 = step(hRSDec,decDataInto2);
        deDataRS3 = step(hRSDec,decDataInto3);

        
        % Convert to the bit form
        deDataRSB1 = de2bi(deDataRS1,8);
        deDataRSB2 = de2bi(deDataRS2,8);
        deDataRSB3 = de2bi(deDataRS3,8);
        
        % Compute the errors for each data with the original data
        BERVec1(:,n) = step(hErrorCalc1,data(:),deDataRSB1(:));
        BERVec2(:,n) = step(hErrorCalc2,data(:),deDataRSB2(:));
        BERVec3(:,n) = step(hErrorCalc3,data(:),deDataRSB3(:));
        BERVec4(:,n) = step(hErrorCalc4,data(:),DemodData4);
    end
end
% save the data
BER_QPSK1 = BERVec1(1,:);
BER_QPSK2 = BERVec2(1,:);
BER_QPSK3 = BERVec3(1,:);
BER_QPSK4 = BERVec4(1,:);
save('QPSK_sim2.mat','EbNoEncoderInput','BER_QPSK1','BER_QPSK2','BER_QPSK3','BER_QPSK4');

%% Plot Results
load('QPSK_sim2.mat');
semilogy(EbNoEncoderInput,BER_QPSK1);
hold on
semilogy(EbNoEncoderInput,BER_QPSK2);
hold on
semilogy(EbNoEncoderInput,BER_QPSK3);
hold on
semilogy(EbNoEncoderInput,BER_QPSK4)
grid on;
title('Example of Communication System Simulator');
ylabel('BER');
xlabel('EbNo Ratio');
legend('1/2','2/3','5/6','uncoded')