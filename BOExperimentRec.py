import pandas as pd
import torch
import numpy as np
from scipy.optimize import minimize
from botorch.models.fully_bayesian import SaasFullyBayesianSingleTaskGP
from botorch.models import ModelListGP
from botorch.models.transforms.input import InputTransform
from botorch.models.transforms.outcome import Standardize
from botorch.fit import fit_fully_bayesian_model_nuts
from torch.distributions import Normal
import warnings

# Suppress warnings
warnings.filterwarnings('ignore')

# Input transform function for GP use
class PhysicsBridge(InputTransform, torch.nn.Module):
    def __init__(self, bounds):
        super().__init__()
        self.transform_on_train = True
        self.xL = np.hstack([[bounds[i][0] for i in range(8)],[5,5,5,5]])
        self.xU = np.hstack([[bounds[i][1] for i in range(8)],[180,180,180,180]])
        
    def forward(self, X):
        # Handles raw 8-dim input and adds the 4 products
        m1a, m1b, m2a, m2b = X[..., 0], X[..., 1], X[..., 2], X[..., 3]
        l1a, l1b, l2a, l2b = X[..., 4], X[..., 5], X[..., 6], X[..., 7]
        p = torch.stack([m1a*l1a, m1b*l1b, m2a*l2a, m2b*l2b], dim=-1)
        Xnorm = (torch.cat([X, p], dim=-1) - self.xL)/(self.xU-self.xL) 
        return Xnorm
    
    def transform(self, X): 
        return self.forward(X)
    
    def untransform(self,X):
        X_raw = X * (self.xU - self.xL) + self.xL
        return X_raw[..., :8]

# Load training data for GP    
def loadData():
    df = pd.read_csv('training_data.csv', skip_blank_lines=True)
    df = df.dropna(how='all')
    X_raw = torch.tensor(df[['M1A','M1B','M2A','M2B','L1A','L1B','L2A','L2B']].values, dtype=torch.double)
    Y = torch.tensor(df[['Y1', 'Y2']].values, dtype=torch.double)
    Y_var = torch.tensor(df[['Y1_std', 'Y2_std']].values, dtype=torch.double)**2

    return X_raw, Y, Y_var

# Train SAASGP
def trainGP(X,Y,Yvar,bridge):
    # Moniter Y1: interface slope    
    gp_Y1 = SaasFullyBayesianSingleTaskGP(train_X=X, 
                                          train_Y=Y[:, [0]], 
                                          train_Yvar=Yvar[:, [0]], 
                                          input_transform=bridge, 
                                          outcome_transform=Standardize(m=1))
    
    # Moniter Y2: registery number    
    gp_Y2 = SaasFullyBayesianSingleTaskGP(train_X=X, 
                                          train_Y=Y[:, [1]], 
                                          train_Yvar=Yvar[:, [1]], 
                                          input_transform=bridge, 
                                          outcome_transform=Standardize(m=1))
    
    fit_fully_bayesian_model_nuts(gp_Y1, warmup_steps=256, num_samples=128, thinning=16, disable_progbar=True)
    fit_fully_bayesian_model_nuts(gp_Y2, warmup_steps=256, num_samples=128, thinning=16, disable_progbar=True)
    model = ModelListGP(gp_Y1, gp_Y2).eval()

    return model

# Acquisition function
def acqfunc(x_numpy,model,best_f):
    '''
    Takes numpy array from Scipy, returns negative CEI (to minimize).
    '''
    # Convert numpy to torch tensor
    xt = torch.from_numpy(x_numpy).to(dtype=torch.double).view(1, 1, 8)
    
    with torch.no_grad():
        # Get predictions (handles the 8->12 transform internally)
        post = model.posterior(xt)
        mu = post.mean # (128, 1, 1, 2)
        sigma = post.variance.clamp_min(1e-9).sqrt()
        
        # Calculate EI for Y1
        mu1, sig1 = mu[..., 0], sigma[..., 0]
        u = (mu1 - best_f) / sig1
        dist = Normal(0, 1)
        ei = (mu1 - best_f) * dist.cdf(u) + sig1 * torch.exp(dist.log_prob(u))
        
        # Calculate Prob(Y2 >= 0.7)
        mu2, sig2 = mu[..., 1], sigma[..., 1]
        pf = dist.cdf((mu2 - 0.7) / sig2)
        
        # Average the 128 MCMC samples into one single score
        cei = (ei * pf).mean().item()
        
    return -cei # Negative because Scipy minimizes

# Get recommendation for the next experiment
def getRec(bnds,model,best_f):
    # Constraints: fun(x) >= 0
    xL = np.array([bnds[i][0] for i in range(8)])
    xU = np.array([bnds[i][1] for i in range(8)])

     # input constraints   
    cons = [
        {'type': 'ineq', 'fun': lambda x: x[0]*x[4] - x[1]*x[5]}, # m1a*l1a >= m1b*l1b
        {'type': 'ineq', 'fun': lambda x: x[3]*x[7] - x[2]*x[6]}  # m2b*l2b >= m2a*l2a
    ]
    
    best_score = float('inf')
    best_x = None
    
    for i in range(25):
        # Random valid start point
        x0 = np.random.rand(8)
        
        res = minimize(
            acqfunc, 
            x0=x0, 
            bounds=[(0,1)]*8, 
            constraints=cons, 
            method='SLSQP',
            options={'maxiter': 100, 'ftol': 1e-4},
            args = (model, best_f)
        )
        
        if res.success and res.fun < best_score:
            best_score = res.fun
            best_x = res.x
        #print(f"Restart {i+1} complete.", end="\r")
    intx = finalize(best_x*(xU-xL)+xL,bnds,model,best_f)
    return intx

# Output function
def finalize(x,bnds,model,best_f):
    '''
    Take the recommendation point from the acquisition func (continuous space) as a base,
    search multiple neighborhood points that satisfied the input constraints (integer space),
    return the candidate list.
    '''
    mdfr = 0.20
    xvalid = False

    xL = np.array([bnds[i][0] for i in range(8)])
    xU = np.array([bnds[i][1] for i in range(8)])

    while xvalid is False:
        # Generate random points within +/- 20% of initial point
        low = np.floor(x-mdfr*xL)
        low = np.clip(low,xL,None)
        high = np.ceil(x+mdfr*xU)
        high = np.clip(high,None,xU)
        npoints = 1000
        
        xgen = np.random.randint(low, high+1, size=(npoints,8))
        xgen = np.clip(xgen,low,high)

        # Ensure physics constraints hold after rounding
        mask1 = np.where((xgen[:,1]*xgen[:,5]) <= (xgen[:,0]*xgen[:,4]))
        mask2 = np.where((xgen[:,2]*xgen[:,6]) <= (xgen[:,3]*xgen[:,7]))
        mask = np.intersect1d(mask1,mask2)

        if len(mask) > 0:
            xvalid = True
        else:
            mdfr = mdfr + 0.05
    
    xpass = xgen[mask,:]
    xscl = (xpass-xL)/(xU-xL)

    bestaf = -np.inf
    bestx = np.zeros([1,8])
    for i in range(len(xscl)):
        xaf = -acqfunc(xscl[i,:],model,best_f)

        if xaf > bestaf:
            bestaf = xaf
            bestx = xscl[i,:]
    
    # Return point closest to rounded point
    #dist = np.linalg.norm(c-xpass,axis=1)
    #bestx = xpass[np.argmin(dist),:].reshape(1,-1)
    return bestx * (xU-xL)+xL

# Generate the recommendation.txt file 
def printNewPoint():
    '''
    Rank the acquisition function value using the candidate list from 'finalize' func,
    return the best candidate.
    '''
    # BO Loop for automatic experiment recommendation
    # Bounds: [M-range, L-range]
    bnds = [(1,3)]*4 + [(5,60)]*4
    bridge = PhysicsBridge(bnds)
    
    # Load complete training data list
    X, Y, Yvar = loadData()
    
    # Find current best Y1 where Y2 >= 0.7
    mask = Y[:, 1] >= 0.7
    best_f = Y[mask, 0].max() if mask.any() else Y[:, 0].min()
    
    # Train GP models
    model = trainGP(X,Y,Yvar,bridge)
    
    # Find optimal experiment
    xrec = getRec(bnds,model,best_f)
    xrec = xrec.reshape(1,-1)

    # Print xrec to txt file [M1A, M1B, M2A, M2B, L1A, L1B, L2A, L2B]
    np.savetxt('recommendation.txt', xrec, fmt = '%3.0i', delimiter=',')
    
if __name__ == "__main__":
    printNewPoint()
