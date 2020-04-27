# Standard imports
import torch 
import numpy as np
import scipy.io   

# Loading the confusion matrix

confusionMatrix = scipy.io.loadmat("../data/brailletouch.mat")
confusionMatrix = confusionMatrix["A"]

# Loading the feature vector matrix according to either of the peano curves 

rowFeatureVectors =  scipy.io.loadmat("../data/FVrow.mat")
rowFeatureVectors = rowFeatureVectors['Xrow']

# Initializing the distance metric 

M = torch.rand(6,6)

print(M)
