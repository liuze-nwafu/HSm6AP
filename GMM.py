def GMM(path_input,file_input,path_output,file_output):
    import numpy as np
    import pandas as pd
    from sklearn.mixture import GaussianMixture
    file = path_input+file_input
    X = pd.read_csv(open(file, "r"), delimiter=",", skiprows=0)
    gmm = GaussianMixture(n_components=5)
    gmm.fit(X)
    labels = gmm.predict(X)
    np.savetxt(path_output+file_output, labels, delimiter=",")
    