function dpcmm()
clc;
clear all;
close all;
  Inp = (imread('lena512.bmp'));
  I = double(Inp);
%  I =double( [92 93 44 98 80;96 54 56 43 41; 97 93 100 45 92; 104 92 83 80 100;56 89 75 68 103])
I = I - 128;
 
[N,M] = size(I);
 
a=1;b=-1;c=1;
  XInputImageInfo = imfinfo('Lena512.bmp');
      
       imwrite(Inp,'JInputImage.jpg');
     JPEGInputImageInfo = (imfinfo('JInputImage.jpg'));
     JPEGImageInfoFileSize = JPEGInputImageInfo.FileSize * ones(1,12)
         JI = imread('JInputImage.jpg');
     XInputImageInfoFileSize = XInputImageInfo.FileSize * ones(1,12)
     psnrr = [];
     snrr = [];
      jpsnrr = [];
     jsnrr = [];
     deltas = [];
     sizee = [];
     qpe = [];
     

for delta = 1:12
 
 E(1,1:M) = I(1,1:M);
    E(2:N,1) = I(2:N,1);
    CR = E;
for i=2:N 
    for j = 2:M
        p = a*CR(i,j-1) + b*CR(i-1,j-1) + c*CR(i-1,j);
        E(i,j) = Q(I(i,j) - p,delta);
        CR(i,j) = Qd(E(i,j),delta) + p;
    end
end 
[ B dict] = entropy_code(E);
lengthOfB = numel(de2bi(B))
 
dict = 0;
 
R = dpcm_decode(E,delta,N,M,a,b,c,dict);
imwrite(R,'barbaraReconstructed.jpg');
    ReconstructedImageInfo = imfinfo('barbaraReconstructed.jpg');
 
     if delta == 12 || delta == 1
        figure('Name','Barbara Image')
        
     imshow([Inp R]);
         title('Original Image - - - - - - - - -  - - - - - - - - - - - - - - - - - - - - - - - - - Reconstructed Image');    

    end

[psnrval, snr] = psnr(R,Inp);
[jpsnrval, jsnr] = psnr(JI,Inp);
    psnrr = [ psnrr psnrval];
    jpsnrr = [jpsnrr jpsnrval];
    jsnrr = [jsnrr jsnr];
    snrr = [snrr snr];
     deltas = [deltas delta];
    sizee = [sizee ReconstructedImageInfo.FileSize];
    qpe = [qpe lengthOfB/8];
end

 
    
    figure()
    plot(deltas,psnrr,'g--o',deltas,snrr,'r:*',deltas,jpsnrr,'m-',deltas,jsnrr,'c-')
    title('Lena Quantization Interval vs SNR and PSNR')
    legend('PSNR of Nearlossless EntropyDPCM','SNR of NearlosslessEntropyDPCM','PSNR of JPEG Image(No quantization)','SNR of JPEG Image(No quantization)');
    ylabel('PSNR and SNR in decibel');
    xlabel('Quantization Interval');
    
    figure()
     plot(deltas,qpe/1000,'r--*',deltas,XInputImageInfoFileSize/1000,'b:',deltas,JPEGImageInfoFileSize/1000,'m-.')
    title('Lena Decoded Image Size using Entropy coded DPCM')
    legend('EntropyDPCM encoded Size','JPEG encoded size','Original bitmap image Size');
    ylabel('Encoded size in Kbytes');
    xlabel('Quantization Interval');
    
     figure()
     plot(deltas,XInputImageInfoFileSize./qpe,'c--o',deltas,XInputImageInfoFileSize./JPEGImageInfoFileSize,'m-')
    title('Lena compression ratio comparison between JPEG and Entropy coded DPCM')
    legend('EntropyDPCM Compression Ratio','JPEG Compression Ratio');
    ylabel('Compression Ratio');
    xlabel('Quantization Interval');
    
    JpegCompressionRatio = XInputImageInfoFileSize./JPEGImageInfoFileSize
    EntropyDPCM = XInputImageInfoFileSize./qpe
end

function quan = Q(I,delta)
quan = round(I/delta) ;
end

function quant = Qd(I,delta)
quant = I * delta;
end

function [B dict] = entropy_code(img)
[M N] = size(img);
numPix = M*N;
symbols = unique(img);
co = 1;
for i =symbols'
     prob(co)=numel(find(img==i));
     co = co + 1;
end
prob = prob/numPix;
[dict,avglen] = huffmandict(symbols,prob); 
 B = huffmanenco(img(:),dict);
end

function D = entropy_decode(img,dict)
dsig = huffmandeco(img,dict);
D = reshape(dsig,N,M);
end

function D = dpcm_decode(img,delta,N,M,a,b,c,dict)
% E = entropy_decode(img,dict);
E = img;
 Inew(1,1:M) = E(1,1:M);
    Inew(2:N,1) = E(2:N,1); 
for i = 2:N
    for j = 2:M
        Inew(i,j) = Qd(E(i,j),delta) + (a*Inew(i,j-1)) + (b*Inew(i-1,j-1)) + (c*Inew(i-1,j));
    end
end
val = Inew + 128;
D = uint8(val);
 
end