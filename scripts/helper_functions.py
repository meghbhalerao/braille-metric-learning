import torch 
import torch.nn
import numpy as np

# Z scoring the matrix 
def normalize(matrix):
    mean = torch.mean(matrix.flatten())
    sigma = torch.std(matrix.flatten())
    matrix = (matrix - mean)/sigma
    return matrix 

# Making a non-symmetrix square matrix symmetric 
def makeSymmetric(matrix):
    return 0.5*(matrix + matrix.T)

# Defining some loss functions here

def MSELoss(mat1,mat2):
    loss = torch.nn.MSELoss()
    return loss(mat1,mat2)
    
                                        
