import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
#import seaborn as sns
import sklearn as sk
import tensorflow as tf
from tensorflow import keras
if not tf.config.list_physical_devices('GPU'):
    print("No GPU was detected. LSTMs and CNNs can be very slow without a GPU.")
import LD_Lib
from sklearn.model_selection import train_test_split

############################### UNSUPERVISED AND DIM REDUCTION ################################

def TSNE(NT, c = None, PLOT_IT = False):
    from sklearn.manifold import TSNE
    tsne = TSNE()
    X_valid_2D = tsne.fit_transform(NT)
    X_valid_2D = (X_valid_2D - X_valid_2D.min()) / (X_valid_2D.max() - X_valid_2D.min())
    if PLOT_IT:
        fig = plt.figure()
        if c is None:
            plt.scatter(X_valid_2D[:, 0], X_valid_2D[:, 1], s=10, cmap="tab10")
        else:
            plt.scatter(X_valid_2D[:, 0], X_valid_2D[:, 1], c=c, s=10, cmap="tab10")

        plt.axis("off")
        #plt.show(block=False)
        #plt.pause(0.1)

def PCA(NT, PLOT_IT = False):
    from sklearn.decomposition import PCA

    pca = PCA()
    PCscores = pca.fit_transform(NT)
    pca.fit(NT)
    cumsum = np.cumsum(pca.explained_variance_ratio_)
    d90 = np.argmax(cumsum >= 0.90) + 1

    if PLOT_IT:
        fig = plt.figure(figsize=(4,3))
        plt.plot(PCscores[:,0], PCscores[:, 1], "b.")
        plt.xlabel("$z_1$", fontsize=18)
        plt.ylabel("$z_2$", fontsize=18, rotation=0)
        plt.grid(True)

        fig = plt.figure()
        plt.plot(PCscores[:,0], "b.")
        plt.plot(PCscores[:,1], "r.")
        plt.xlabel("time", fontsize=18)
        plt.ylabel("pca", fontsize=18, rotation=0)
        plt.grid(True)

        plt.figure(figsize=(6,4))
        plt.plot(cumsum, linewidth=3)
        plt.axis([0, 400, 0, 1])
        plt.xlabel("Dimensions")
        plt.ylabel("Explained Variance")
        plt.plot([d90, d90], [0, 0.90], "k:")
        plt.plot([0, d90], [0.90, 0.90], "k:")
        plt.plot(d90, 0.90, "ko")
        #plt.annotate("Elbow", xy=(65, 0.85), xytext=(70, 0.7),
        #            arrowprops=dict(arrowstyle="->"), fontsize=16)
        plt.grid(True)

        # plt.show()

    return(PCscores, pca, d90)

def Kernel_PCA(NT, PLOT_IT = False):
    # incredibly slow
    from sklearn.decomposition import KernelPCA

    lin_pca = KernelPCA(n_components = 2, kernel="linear", fit_inverse_transform=True)
    rbf_pca = KernelPCA(n_components = 2, kernel="rbf", gamma=0.0433, fit_inverse_transform=True)
    sig_pca = KernelPCA(n_components = 2, kernel="sigmoid", gamma=0.001, coef0=1, fit_inverse_transform=True)

    t = range(0,NT.shape[0])
    
    plt.figure(figsize=(11, 4))
    for subplot, pca, title in ((131, lin_pca, "Linear kernel"), (132, rbf_pca, "RBF kernel, $\gamma=0.04$"), (133, sig_pca, "Sigmoid kernel, $\gamma=10^{-3}, r=1$")):
        X_reduced = pca.fit_transform(NT)
        if subplot == 132:
            X_reduced_rbf = X_reduced
        
        plt.subplot(subplot)
        #plt.plot(X_reduced[y, 0], X_reduced[y, 1], "gs")
        #plt.plot(X_reduced[~y, 0], X_reduced[~y, 1], "y^")
        plt.title(title, fontsize=14)
        plt.scatter(X_reduced[:, 0], X_reduced[:, 1], c=t, cmap=plt.cm.hot)
        plt.xlabel("$z_1$", fontsize=18)
        if subplot == 131:
            plt.ylabel("$z_2$", fontsize=18, rotation=0)
        plt.grid(True)

    # plt.show()

def Kmeans(NT, k_range = range(2,10), PLOT_IT = False):
    from sklearn.cluster import KMeans
    kmeans_per_k = [KMeans(n_clusters=k, random_state=42).fit(NT) for k in k_range]
    inertias = [model.inertia_ for model in kmeans_per_k]
    if PLOT_IT:
        fig = plt.figure(figsize=(4,3))
        plt.plot(inertias, "b.")
        plt.grid(True)
        #plt.show()

    return(kmeans_per_k)

def Autoencode_1Layer(NT, PLOT_IT = False):
    # Python ≥3.5 is required
    # Scikit-Learn ≥0.20 is required
    # Build the network.
    
    n_time_bins = NT.shape[0]
    n_neurons = NT.shape[1]
    RIX = np.arange(0, n_time_bins)    # array of all indices
    np.random.shuffle(RIX) 
    mid = np.round(n_time_bins/2).astype(np.int16)
    X_train, X_valid = NT[RIX[:-mid],:], NT[RIX[-mid:],:]
    X_train = X_train - X_train.mean(axis=0, keepdims=0)
    
    # Create the network.
    encoder = keras.models.Sequential([keras.layers.Dense(2, input_shape=[n_neurons])])
    decoder = keras.models.Sequential([keras.layers.Dense(n_neurons, input_shape=[2])])
    autoencoder = keras.models.Sequential([encoder, decoder])

    autoencoder.compile(loss="mse", optimizer=keras.optimizers.SGD(lr=1.5))
    history = autoencoder.fit(X_train, X_train, epochs=20)
    codings = encoder.predict(X_train)
    if PLOT_IT:
        fig = plt.figure(figsize=(4,3))
        plt.plot(codings[:,0], codings[:, 1], "b.")
        plt.xlabel("$z_1$", fontsize=18)
        plt.ylabel("$z_2$", fontsize=18, rotation=0)
        plt.grid(True)
        #plt.show(block=False)

    return(codings,encoder)

def Autoencode_Stacked(NT, n_units = (8, 8, 6), epochs = 50, PLOT_IT = False):
    # Python ≥3.5 is required
    # Scikit-Learn ≥0.20 is required
    # Build the network.
    #from keras.regularizers import l1
    from tensorflow.keras.regularizers import l1
    n_time_bins = NT.shape[0]
    n_neurons   = NT.shape[1]
    RIX = np.arange(0, n_time_bins)    # array of all indices
    np.random.shuffle(RIX) 
    mid = np.round(n_time_bins/2).astype(np.int16)
    X_train, X_valid = NT[RIX[:-mid],:], NT[RIX[-mid:],:]
    X_train = X_train - X_train.mean(axis=0, keepdims=0)
    #LD_Lib.Plot_Matrix(X_train[0:10000,:])
    #plt.show()
    # Create the network.
    #tf.random.set_seed(42)
    #np.random.seed(42)

    #def rounded_accuracy(y_true, y_pred):
    #    return keras.metrics.binary_accuracy(tf.round(y_true), tf.round(y_pred))
    
    encoder = keras.models.Sequential([
        keras.layers.Dense(n_units[0], input_shape=[n_neurons]),
        keras.layers.Dense(n_units[1], activation="selu"), # , activity_regularizer=l1(0.001)
        keras.layers.Dense(n_units[2], activation="selu"),
    ])
    # NOTE: output needs to be linear for spikes and autoencoder and mse seemed to work best and not the rounded_accuracy.
    decoder = keras.models.Sequential([
        keras.layers.Dense(n_units[1], activation="selu", input_shape=[n_units[2]]),
        keras.layers.Dense(n_units[0], activation="selu"),
        keras.layers.Dense(n_neurons, activation="linear"),
    ])
    stacked_ae = keras.models.Sequential([encoder, decoder]) 
    #stacked_ae.compile(loss="binary_crossentropy",
    #                optimizer=keras.optimizers.SGD(lr=1.5), metrics=[rounded_accuracy])
    stacked_ae.compile(loss="mse")
    history = stacked_ae.fit(X_train, X_train, epochs=epochs, validation_data=(X_valid, X_valid))

    codings_orig = encoder.predict(NT)

    if PLOT_IT:
        codings = encoder.predict(X_train)
        codings_valid = encoder.predict(X_valid)
        Plot_Learning_Curves(history.history["loss"], history.history["val_loss"])
        Plot_Layers(encoder.layers)
        pd.plotting.scatter_matrix(pd.DataFrame(codings),hist_kwds={'bins': 40})
        #plt.figure(202)
        #sns.pairplot(pd.DataFrame(codings),markers=".",kind="reg",diag_kind="kde")
    #    return(encoder, codings_orig, history)

    return(codings_orig, encoder,stacked_ae)

############################### Classificiation and regression ################################
def Linear_Regression(NT,Y,PLOT_IT = False):
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import mean_squared_error

    n_time_bins = NT.shape[0]
    n_neurons   = NT.shape[1]
    RIX = np.arange(0, n_time_bins)    # array of all indices
    np.random.shuffle(RIX) 
    mid = np.round(n_time_bins/2).astype(np.int16)
    X_train, X_valid = NT[RIX[:-mid],:], NT[RIX[-mid:],:]
    #X_train = X_train - X_train.mean(axis=0, keepdims=0)
    Y_train, Y_valid = Y[RIX[:-mid],:], Y[RIX[-mid:],:]

    reg = LinearRegression()
    reg.fit(X_train, Y_train)
    Y_pred = reg.predict(X_valid)
    Y_pred[Y_pred>1000] = 0
    Y_pred[Y_pred<0] = 0
    mse = mean_squared_error(Y_valid, Y_pred)
    if PLOT_IT:
        Plot_Pred_Vs_Valid(Y_pred, Y_valid,'Lin Regress mse ' + f"{mse:.1f}")

    return(mse, reg)


def Decision_Tree_Regression(NT,Y,PLOT_IT = False):
    from sklearn.metrics import mean_squared_error
    from sklearn.tree import DecisionTreeRegressor

    import scipy

    n_time_bins = NT.shape[0]
    n_neurons   = NT.shape[1]
    RIX = np.arange(0, n_time_bins)    # array of all indices
    np.random.shuffle(RIX) 
    mid = np.round(n_time_bins/2).astype(np.int16)
    X_train, X_valid = NT[RIX[:-mid],:], NT[RIX[-mid:],:]
    #X_train = X_train - X_train.mean(axis=0, keepdims=0)
    Y_train, Y_valid = Y[RIX[:-mid]], Y[RIX[-mid:]]

    reg = DecisionTreeRegressor(random_state=42)
    reg.fit(X_train, Y_train)
    Y_pred = reg.predict(X_valid)
    Y_pred[Y_pred>1000] = 0
    Y_pred[Y_pred<0] = 0
    mse = mean_squared_error(Y_valid, Y_pred)
    if PLOT_IT:
        Plot_Pred_Vs_Valid(Y_pred, Y_valid,'Tree Regress mse ' + f"{mse:.1f}")

    return(mse, reg)

def SVM(NT, Y, kernel = 'rbf', PLOT_IT = False):
    from sklearn.svm import SVR
    from sklearn.metrics import mean_squared_error
    import scipy

    n_time_bins = NT.shape[0]
    n_neurons   = NT.shape[1]
    RIX = np.arange(0, n_time_bins)    # array of all indices
    np.random.shuffle(RIX) 
    mid = np.round(n_time_bins/2).astype(np.int16)
    X_train, X_valid = NT[RIX[:-mid],:], NT[RIX[-mid:],:]
    #X_train = X_train - X_train.mean(axis=0, keepdims=0)
    Y_train, Y_valid = Y[RIX[:-mid]], Y[RIX[-mid:]]

    reg = SVR(kernel=kernel)
    reg.fit(X_train, Y_train)
    Y_pred = reg.predict(X_valid)
    mse = mean_squared_error(Y_valid, Y_pred)
    rmse = np.sqrt(mse)

    if PLOT_IT:
        Plot_Pred_Vs_Valid(Y_pred, Y_valid,'SVM mse ' + f"{mse:.1f}")


    return(rmse, reg)

def MLP(NT, Y , PLOT_IT = False):
    from sklearn.metrics import mean_squared_error
    import scipy

    IXtrall, IXtest, IXtr, IXval = Split_Train_Test(NT, 0.2, 0.3)
    X_train,Y_train = NT[IXtr,:],Y[IXtr]
    X_valid,Y_valid = NT[IXval,:],Y[IXval]
    X_test,Y_test = NT[IXtest,:],Y[IXtest]
    
    model = keras.Sequential()
    model.add(keras.layers.Dense(20, input_shape=(X_train.shape[1],)))
    model.add(keras.layers.Dense(1))

    # model = keras.models.Sequential([
    #     keras.layers.Dense(30, activation="relu", input_shape=X_train.shape[1]),
    #     keras.layers.Dense(1)
    # ])

    model.compile(loss="mean_squared_error", optimizer=keras.optimizers.SGD(lr=1e-3))
    #model.fit(x, y, batch_size=32, epochs=10)
    history = model.fit(X_train, Y_train, epochs=20, validation_data=(X_valid, Y_valid))

    Y_pred = model.predict(X_test)

    mse = model.evaluate(X_test, Y_test)

    if PLOT_IT:
        Plot_Pred_Vs_Valid(Y_pred, Y_test,'MLP mse(test) ' + f"{mse:.1f}")

    return(mse, model)


def LSTM(NT, Y , PLOT_IT = False):
    # This is not working yet as it does not work with a continuous outcome vbl. I think i will have to discritize the locatio vbl.
    from sklearn.metrics import mean_squared_error
    import scipy
    from sklearn.naive_bayes import GaussianNB

    IXtr, IXtest,a,b = Split_Train_Test(NT, test_ratio = 0.4)
    X_train,Y_train = NT[IXtr,:],Y[IXtr]
    X_test,Y_test = NT[IXtest,:],Y[IXtest]
    

    model = keras.models.Sequential([
        keras.layers.LSTM(20, return_sequences=True, input_shape=[None, 1]),
        keras.layers.LSTM(20, return_sequences=True),
        keras.layers.TimeDistributed(keras.layers.Dense(10))
    ])

    model.compile(loss="mse", optimizer="adam", metrics=[last_time_step_mse])
    history = model.fit(X_train, Y_train, epochs=20,
                        validation_data=(X_valid, Y_valid))
                        
    history = model.fit(X_train, Y_train.squeeze())
    Y_pred = model.predict(X_test)

    mse = model.evaluate(X_test, Y_test)

    if PLOT_IT:
        plt.figure()
        Plot_Learning_Curves(history.history["loss"], history.history["val_loss"])
  
        Plot_Pred_Vs_Valid(Y_pred, Y_test,'MLP mse ' + f"{mse:.1f}")

    return(mse, model)


def Naive_Bayes_Gaussian(NT, Y , PLOT_IT = False):
    # This is not working yet as it does not work with a continuous outcome vbl. I think i will have to discritize the locatio vbl.
    from sklearn.metrics import mean_squared_error
    import scipy
    from sklearn.naive_bayes import GaussianNB

    IXtr, IXtest,a,b = Split_Train_Test(NT, test_ratio = 0.4)
    X_train,Y_train = NT[IXtr,:],Y[IXtr]
    X_test,Y_test = NT[IXtest,:],Y[IXtest]
    
    model = GaussianNB()
    #model.fit(x, y, batch_size=32, epochs=10)
    history = model.fit(X_train, Y_train.squeeze())
    Y_pred = model.predict(X_test)

    mse = model.evaluate(X_test, Y_test)

    if PLOT_IT:
        plt.figure()
        Plot_Learning_Curves(history.history["loss"], history.history["val_loss"])

        Plot_Pred_Vs_Valid(Y_pred, Y_test,'MLP mse ' + f"{mse:.1f}")

    return(mse, model)

############################# GENERAL UTILITIES #####################################
def Validate_Against_Behavior(SC,COMP, PLOT_IT = False):
    # Assume that SC is, for example, some new dimensions that describe the data and COMP are data that you wish to compare.
    # Compute measures of similarity for each combination of SC to COMP.
    # The easiest metric is r-squared.
    import scipy.stats
    r = np.zeros((SC.shape[1],COMP.shape[1]))
    p = np.zeros((SC.shape[1],COMP.shape[1]))

    for iSC in np.arange(0,SC.shape[1]):
        for iCOMP in np.arange(0,COMP.shape[1]):
            rp = scipy.stats.pearsonr(SC[:,iSC], COMP[:,iCOMP]) # r[iSC,iCOMP] = 
            r[iSC,iCOMP] = rp[0]
            p[iSC,iCOMP] = rp[1]
    if PLOT_IT:
        plt.figure()
        plt.subplot(1,2,1)
        plt.pcolor(r)
        plt.colorbar()
        plt.title('r')
        plt.subplot(1,2,2)
        plt.pcolor(r*p<0.05)
        plt.title('p<0.05')
        #plt.colorbar()

    return(r,p)

def Split_Train_Test(NT, test_ratio = 0.5, test_valid_ratio = 0.5):
    shuffled_indices = np.random.permutation(NT.shape[0])
    test_set_size = int(NT.shape[0] * test_ratio)
    test_indices = shuffled_indices[:test_set_size]
    all_train_indices = shuffled_indices[test_set_size:]
    # Further split the train into a train and validation set.
    shuffled_indices = np.random.permutation(len(all_train_indices))
    valid_set_size = int(len(all_train_indices) * test_valid_ratio)
    valid_indices = shuffled_indices[:valid_set_size]
    test_indices = shuffled_indices[valid_set_size:]

    return all_train_indices, test_indices, test_indices, valid_indices


############################# PLOTTING #####################################

def Plot_Learning_Curves(loss, val_loss):
    # Call as plot_learning_curves(history.history["loss"], history.history["val_loss"])
    plt.figure()
    plt.plot(np.arange(len(loss)) + 0.5, loss, "b.-", label="Training loss")
    plt.plot(np.arange(len(val_loss)) + 1, val_loss, "r.-", label="Validation loss")
    #plt.gca().xaxis.set_major_locator(plt.ticker.MaxNLocator(integer=True))
    plt.axis([1, 20, 0, 0.05])
    plt.legend(fontsize=14)
    plt.xlabel("Epochs")
    plt.ylabel("Loss")
    plt.grid(True)

def Plot_Layers(layers):
    # expects model.layers input
    plt.figure()

    cnt = 0
    for layer in layers:
        wt = layer.get_weights() # list of numpy arrays
        for thewt in wt:
            cnt = cnt + 1
            plt.subplot(len(layers),2,cnt)
            if thewt.ndim == 1:
                plt.plot(thewt)
            else:
                plt.pcolor(thewt)
                plt.colorbar()

def Plot_Pred_Vs_Valid(Y_pred, Y_valid, title_str):
    import scipy
    #from sklearn.metrics import mean_squared_error
    plt.figure()

    r = scipy.stats.linregress(Y_pred.squeeze(), Y_valid.squeeze())
    plt.subplot(1,2,1)
    plt.plot(Y_pred,Y_valid,'.')
    plt.plot(Y_pred, r.intercept + r.slope*Y_pred, 'r', label='fitted line')
    plt.title(title_str)

    plt.subplot(1,2,2)
    plt.plot(Y_pred,'.')
    plt.plot(Y_valid,'.')
    plt.title(title_str)