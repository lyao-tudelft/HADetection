# HADetection
Human activity detection project for Wireless Networking.

File “main.m” is the master program simulating a mutipath channel attenuation, while other files are functions fullfilling specific functionalities serving for the master program. To see how the program goes, run "main.m".

Folder "SAGE.m" contains codes calculating channel state information (CSI) out of multipath attenuated signals. File "test" generates trasmitted signal , calculates the signal and multipath corrupted output signal from the channel, then does parameter estimate to the channel by SAGE and Wiener Filter. Function "SAGEinit.m" does initialization of estimate, function "SAGE.m" does main phase of estimate. File "correlation.m" does post processing to estimated channel CSI to extract movement out of the time-varying channel.
