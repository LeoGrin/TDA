import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def PCA(P, normalize, k, names):
    # Preprocessing
    P = P - np.mean(P, axis=0)
    if normalize:
        P = P / np.std(P, axis=0)


    n = P.shape[0] # number of points
    d = P.shape[1]

    # Compute the eigenvectors of the covariance matrix
    cov_matrix = (1 / n) * np.dot(P.T, P)
    eigenvalues, eigenvectors = np.linalg.eig([cov_matrix])
    eigenvalues = eigenvalues.reshape(-1)
    indices = np.argsort(eigenvalues)[::-1].reshape(-1)

    #plt.plot([i for i in range(len(indices))], [eigenvalues[i] for i in indices])
    #plt.show()

    # Compute the new points
    proj_matrix = eigenvectors[:, indices[:k]].reshape(d, k)

    new_P = np.dot(P, proj_matrix)
    fig, ax = plt.subplots()
    x = [new_P[i][0] for i in range(n)]
    y = [new_P[i][1] for i in range(n)]
    ax.scatter(x, y)
    for i in range(n):
        ax.annotate(names[i], (x[i], y[i]))

    plt.show()









if __name__ == "__main__":
    df = pd.read_table("decathlon.dat.txt", sep=" ")
    names = df.index.values
    P = df.values
    print(P.shape)
    PCA(P, True, 2, names)
