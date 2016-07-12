# digital_pulse_code_modulation
Implementation of a Near lossless DPCM with analysis to JPEG

In this analysis, A nearloss DPCM compression alone is benchmarked to JPEG in MATLAB (MATLAB implicitly uses Discrete Cosine Transform), 
Note that DPCM is usually used as an intermediary step in compression and can also be used as an intermediary stage in JPEG compression. 
Hence it is not expected to have a greater compression ratio than JPEG.
This analysis just shows the capability of DPCM, The trade off between compression ratio and signal to noise ratio, 
variation as quantization is increased and uses JPEG as a benchmark.
It is also lossy because quantization has been used.
