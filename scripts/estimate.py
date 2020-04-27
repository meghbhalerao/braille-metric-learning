# Standard imports
import torch 
import numpy as np
import scipy.io   


# Z scoring the matrix 
def normalize(matrix):
        mean = torch.mean(matrix.flatten())
        sigma = torch.std(matrix.flatten())
        matrix = (matrix - mean)/sigma
        return matrix 

# Making a non-symmetrix square matrix symmetric 
def makeSymmetric(matrix):
    return 0.5*(matrix + matrix.T)


    

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

res = torch.mm(torch.mm(rowFeatureVectors.T,M.double()),rowFeatureVectors)


for ep in range(num_epochs):
    
    


 
