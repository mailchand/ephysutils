# ephysutils
This toolkit is a set of scripts/functions I developed based on other papers to analyse common electrophysiological indices. There are some other files from other people who have written wonderful tools to perform analysis of electrophysiological data (mainly spike trains) as well as LFP signals.

This is a work in progress that I hope to populate in the next days.

computeLambdaAndBeta(Spks, numIter)

The first function I wrote estimates for a neuron the normalized latency (&#955;) and covariance term (&#946;) parameters described in DiCarlo and Maunsell 2005. The code is somewhat slow but does the job quite well. I implemented the least squares version of this technique recently and found it to be very handy. This may be of utility to some one else if they decide to characterize neural populations in a sensorimotor processing pathway.


