# HADetection
Human activity detection project for Wireless Networking.

File “main.m” is the master program simulating a mutipath channel attenuation, while other files are functions fullfilling specific functionalities serving for the master program. To see how the program goes, run "main.m".

Folder "SAGE.m" contains codes calculating channel state information (CSI) out of multipath attenuated signals. File "test" generates trasmitted signal , calculate the signal and multipath corrupted output signal from the channel, then do parameter estimate to the channel. Function "SAGEinit.m" do initialization of estimate, function "SAGE.m" do main phase of estimate. File "correlation.m" does post processing to estimated channel CSI to extract movement out of the time-varying channel.
