import numpy as np


def modichol(A,alpha,beta):
    n = A.shape[1] # size of A
    L = np.identity(n)
    #################### EXTRA ???
    D = np.zeros((n,1))
    c = np.copy(A) #????
    ######################
    D[0] = np.max(np.abs(A[0,0]),alpha)
    c[:,0]=A[:,0]
    L[1:n,0]=c[1:n,0]/D[0];
    
    for j in range(1,n-1):
        c[j,j]=A[j,j]-(L[j,0:j-1]**2)*D[0:j-1].T
        for i in range(j+1,n):
            c[i,j]=A[i,j]-(   L[i,0:j-1]*L[j,0:j-1]  )*D[0:j-1].T # L[i,0:j-1].*L[j,0:j-1] ??
        theta = np.max(c[j+1:n,j])
        D[j]=np.concatenate(((theta/beta)**2,np.abs(c[j,j]),alpha)).max() 
        L[j+1:n,j]=c[j+1:n,j]/D[j]

    j=n;    
    c[j,j]=A[j,j]-(L[j,0:j-1]**2)*D[0:j-1].T
    D[j]=np.max(np.abs(c[j,j]),alpha) 
    return L*diag(D)*L.T

A = np.array([[1.5756, 1.4679,0.4592], [1.4679, 1.5194, 0.7003], [0.4592, 0.7003,0.9425]]
             
A = np.random.rand(4,4)
A = A*A.T
A1 = modichol(A,.3,0.1)
print(A1)

             
#A = [    1.5756    1.4679    0.4592;
#             1.4679    1.5194    0.7003;
#             0.4592    0.7003    0.9425]
#A1 = modichol(A,0.5,0.5)
# OUTPUT
# A1 = [1.5756    1.4679    0.4592;
#      1.4679    1.8675    0.7003;
#       0.4592    0.7003    0.9425]
             
# function [A1] = modichol(A,alpha,beta )
# %Find a positive definite approximate of symmetric matrix A for Newton method
# Author: Sheragim, 2016 
# n=size(A,2);
# L=eye(n);
# D(1)=max(abs(A(1,1)),alpha);
# c(:,1)=A(:,1);
# L(2:n,1)=c(2:n,1)/D(1);

# for j=2:n-1
#     c(j,j)=A(j,j)-(L(j,1:j-1).^2)*D(1:j-1).';
#     for i=j+1:n
#         c(i,j)=A(i,j)-(L(i,1:j-1).*L(j,1:j-1))*D(1:j-1)';
#     end
#     theta=max(c(j+1:n,j));
#     D(j)=max([(theta/beta)^2,abs(c(j,j)),alpha]);
#     L(j+1:n,j)=c(j+1:n,j)/D(j);
# end
#     j=n;
#     c(j,j)=A(j,j)-(L(j,1:j-1).^2)*D(1:j-1).';
#     D(j)=max(abs(c(j,j)),alpha);
#     A1=L*diag(D)*L';
# end
