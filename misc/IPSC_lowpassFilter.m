% 4-pole butterworth, fc = 1000 Hz
G = [0.0038;
	 0.0037;
	 0.0036;
	 0.0035;
	 1.0000];
SOS = [1.0000    2.0000    1.0000    1.0000   -1.9369    0.9523;
       1.0000    2.0000    1.0000    1.0000   -1.8551    0.8698;
       1.0000    2.0000    1.0000    1.0000   -1.7970    0.8112;
       1.0000    2.0000    1.0000    1.0000   -1.7670    0.7811];