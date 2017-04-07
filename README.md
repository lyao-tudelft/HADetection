# HADetection
Human activity detection project for Wireless Networking.

File “main” is the master program simulating a mutipath channel attenuation, while other files are functions fullfilling specific functionalities serving for the master program. To see how the program goes, run "main".

Folder "SAGE" contains codes calculating channel state information (CSI) out of multipath attenuated signals. File "test" generates trasmitted signal , calculate the signal and multipath corrupted output signal from the channel, then do parameter estimate to the channel. Function "SAGEinit" do initialization of estimate, function "SAGE" do main phase of estimate.
