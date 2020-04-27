# Standard imports
import torch 
import torch.optim as optim
import numpy as np
import scipy.io   
from helper_functions import * 


# Loading the confusion matrix

confusionMatrix = scipy.io.loadmat("../data/brailletouch.mat")
confusionMatrix = torch.tensor(confusionMatrix["A"].astype(float))
confusionMatrix = normalize(confusionMatrix)
confusionMatrix  = makeSymmetric(confusionMatrix)

# Loading the feature vector matrix according to either of the peano curves 

rowFeatureVectors =  scipy.io.loadmat("../data/FVrow.mat")
rowFeatureVectors = torch.tensor(rowFeatureVectors['Xrow'].astype(float))

# Initializing the distance metric 
M = torch.rand(6,6)

# Setting the number of epochs 
num_epochs =  100
M.requires_grad = True


# Setting the ground truth variable as the confusion matrix 
groundTruth = confusionMatrix

# Settin the optimizer to be used to optimize the parameters of the distance matrix
optimizer = optim.Adam([M], lr = 0.1, betas = (0.9,0.999), weight_decay = 0.00005)


for ep in range(num_epochs):
    output = torch.mm(torch.mm(rowFeatureVectors.T,M.double()),rowFeatureVectors)
    loss = MSELoss(output,groundTruth)
    
    print(loss)
    
    
    
    


 
