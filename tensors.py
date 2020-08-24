# -*- coding: utf-8 -*-
"""
Some additional functions for tensor calculus.
"""
import numpy as np


def kelvin(T):
    """
    Transform a 2nd or 4th order tensor to Kelvin notation.

    Parameters
    ----------
    T : numpy.ndarray
        Either a 2nd order tensor (3x3) or a 4th order tensor (3x3x3x3).
        See **Notes** for a discussion of symmetry prerequisites.

    Returns
    -------
    TKel : numpy.ndarray
        Tensor `T` in Kelvin notation. 6x1 if `T` is 3x3, 6x6 if `T` is
        3x3x3x3.
        
    Notes
    -----
    Kelvin notation assumes the tensor :math:`T` exhibits certain symmetry
    properties:
    
    - 2nd order tensors
        2nd order tensors must be symmetric, i.e. :math:`T_{ij}=T_{ji}` .
    
    - 4th order tensors
        4th order tensors must exhibit minor symmetries, i.e.
        :math:`T_{ijkl}=T_{jikl}=T_{ijlk}=T_{jilk}`.
        
        If :math:`T` also exhibits major symmetries, i.e.
        :math:`T_{ijkl}=T_{klij}`, the resulting matrix in Kelvin notation will
        be symmetric.

    """
    TKel = 0
    index = np.array([[0,0],[1,1],[2,2],[0,1],[1,2],[0,2]])
    if T.ndim == 2:
        TKel = np.zeros(6)
        for i in range(6):
           TKel[i] = T[index[i,0], index[i,1]]
           
        TKel[3:6] = TKel[3:6]*np.sqrt(2)

    elif T.ndim == 4:
        TKel = np.zeros((6,6))
        for i in range(6):
            for j in range (6):
                TKel[j,i] = T[index[i,0],index[i,1],index[j,0],index[j,1]]
                
                
        TKel[3:6,3:6] = TKel[3:6,3:6]*np.sqrt(2) 
    
    return TKel


    
def unkelvin(TKel):
    """
    Recover a 2nd or 4th order tensor from Kelvin notation.
    
    Parameters
    ----------
    TKel : numpy.ndarray
        Either a 6x1 vector or a 6x6 matrix, representing a 2nd or 4th order
        tensor, respectively.
        
    Returns
    -------
    T : numpy.ndarray
        Tensor `T` in standard notation. 3x3 is `T` is a 2nd order tensor,
        3x3x3x3 if `T` is a fourth order tensor.
        
    See Also
    --------
    kelvin : For a discussion of Kelvin notation and symmetry prerequisites.
    
    """
    index = np.array([[0,0],[1,1],[2,2],[0,1],[1,2],[0,2]])
    
    if TKel.ndim == 1:
        T = np.zeros((3,3))
        for i in range (6):
            T[index[i,0], index[i,1]] = TKel[i]
            
        
    
    elif TKel.ndim == 2:
        T = np.zeros((3,3,3,3))
        for i in range(6):
           for j in range(6):
               T[index[i,0],index[i,1],index[j,0],index[j,1]] = TKel[j,i]
            
    return T



def identityTensors():
    """
    Compute 4th order identity tensors `IdyI` and `Isym`.

    Returns
    -------
    (IdyI, Isym) : (numpy.ndarray, numpy.ndarray)
        Identity tensors with ``shape=(3,3,3,3)``. See the notes section for
        tensor definitions.
        
    Notes
    -----
    - Index notation for `IdyI`
        :math:`[I \otimes I]_{ijkl} = \delta_{ij} \delta_{kl}`
    - Index notation for `Isym`
        :math:`[I^{sym}]_{ijkl} = \\frac{1}{2} (\delta_{ik} \delta{jl}
        + \delta_{ik} \delta{jl})`
    """
    I = np.eye(3)
    IodyI = np.zeros ((3,3,3,3))
    IudyI = np.zeros ((3,3,3,3))
    Isym = np.zeros ((3,3,3,3))
    IdyI = np.zeros ((3,3,3,3))
    
    for i in range (3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    IodyI[i,j,k,l] = I[i,k]*I[j,l] #????????????????????????
                    Isym[i,j,k,l] = Isym[i,j,k,l]+0.5*(I[i,k]*I[j,l]+I[i,l]*I[j,k]) 
                    IudyI[i,j,k,l] = I[i,l]*I[j,k]
                    IdyI[i,j,k,l] = I[i,j]*I[k,l]
                    
                    
    #Isym =   .5*(IodyI + IudyI)    
    return IdyI, Isym
    
    
def splitVolDev(T):
    """
    Perform a volumetric-deviatoric split of the 2nd order tensor T.
    
    Arguments
    ---------
    T : numpy.ndarray
        If ``T.shape == (3,3)``, it is assumed that `T` is a 2nd order tensor
        in standard notation. If ``T.shape == (6,)``, is is assumend that `T`
        is a symmetric 2nd order tensor in Kelvin notation.
        
    Returns
    -------
    TVol, TDev : (np.ndarray, np.ndarray)
        The volumetric and deviatoric parts of the tensor are returned in the
        same notation as the input, i.e. if `T` was in Kelvin notation, `TVol`
        and `TDev` are in Kelvin notation too.
    """
    pass # Replace this with your own code!

    
#def dyad(a,b):
    """Dyadic product of two 2nd order tensors."""
    #c = ...
                    
    #return c
        