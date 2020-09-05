# Standard imports
import torch 
import torch.optim as optim
import numpy as np
import scipy.io   
from helper_functions import * 


# Loading the confusion matrix

confusionMatrix = scipy.io.loadmat("../data/brailletouch.mat")

# The confusion matrix is converted to a tensor here for further calculations

confusionMatrix = torch.tensor(confusionMatrix["A"].astype(float))
confusionMatrix = normalize(confusionMatrix)

distanceMatrix = 1 - confusionMatrix
distanceMatrix = z_scoring(distanceMatrix)
distanceMatrix  = makeSymmetric(distanceMatrix)

# Loading the feature vector matrix according to either of the peano curves 

rowFeatureVectors =  scipy.io.loadmat("../data/FVrow.mat")
rowFeatureVectors = torch.tensor(rowFeatureVectors['Xrow'].astype(float))

# Initializing the distance metric 
M = torch.rand(6,6)

# Setting the number of epochs 
num_epochs =  10000
M.requires_grad = True


# Setting the ground truth variable as the confusion matrix 
groundTruth = distanceMatrix

# Setting different optimizer to be used to optimize the parameters of the distance matrix
optimizer_adam = optim.Adam([M], lr = 0.0001, betas = (0.9,0.999), weight_decay = 0.00005)
optimizer_sgd = optim.SGD([M], lr = 0.0001)

# Defining variables to store the logging values of during optimizing 
best_loss = 100000
best_epoch = 0

for ep in range(num_epochs):
    output = torch.mm(torch.mm(rowFeatureVectors.T,M.double()),rowFeatureVectors)
    loss = MSELoss(output,groundTruth)
    loss.backward()
    optimizer_sgd.step()
    if(loss.item()<best_loss):
        best_loss = loss.item()
        best_estimate = output
        best_epoch  = ep
    print("The loss for ",ep, "iteration is: ",loss.item())

print("The best loss is at ", best_epoch, "iteration and the best loss is: ", best_loss)
    
print("The matrix which indicates the difference between the distance matrix and estimated distance matrix is:", (best_estimate - distanceMatrix).detach().numpy())
    
    


 
test canges
