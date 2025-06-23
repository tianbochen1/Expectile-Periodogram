import tensorflow as tf
from tensorflow.keras import layers
import numpy as np
from tensorflow.keras import Input
from tensorflow.keras.layers import Lambda
from tensorflow.keras import optimizers
from tensorflow.keras import Model
from tensorflow.keras.utils import to_categorical
from tensorflow.keras.regularizers import l2
import matplotlib.pyplot as plt
from keras import backend as K
import csv
import numpy as np
import time 
from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession
import os
import pandas as pd
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import confusion_matrix
##################LOAD DATA###################################
data = pd.read_csv('ep500.csv')
del data[data.columns[0]]
data = np.array(data)
x = data.reshape( [-1, 46, 500, 1])
y = np.ones([2000,2])
for i in range(1000):
    y[i] = [1,0]
    y[i+1000] = [0,1]
cone = []
hise = []
pree = []
rece = []  
for i in range(10):    
    np.random.seed(i+1)
    z = np.random.permutation(2000) ##shuffle
    x = x[z,]
    x = x.reshape( [-1, 46, 500, 1])
    y = y[z]
    train_ind = np.append(np.arange(800),np.arange(800)+1000)
    test_ind = np.append(np.arange(200)+800,np.arange(200)+1800)
    x_train = x[train_ind,]
    y_train = y[train_ind,]
    x_test = x[test_ind,]
    y_test = y[test_ind,]
    
    
    ##########use GPU#############################################
    #if CPU only, comment the above lines#########################
    os.environ['CUDA_VISIBLE_DEVICES'] = "0"
    config = ConfigProto()
    config.gpu_options.per_process_gpu_memory_fraction = 0.5
    session = InteractiveSession(config=config)
    
    
    #############model############################################
    inp = Input(shape=(46,500,1),dtype = 'float32')
    conv1 = layers.Conv2D(4,(5,5), activation = 'relu',
                        input_shape = (46,500,1),
                        padding = 'same')(inp)
    pool1 = layers.MaxPooling2D((2, 2))(conv1)
    conv2 = layers.Conv2D(8,(5,5), activation = 'relu',
                        padding = 'same')(pool1)
    pool2 = layers.MaxPooling2D((2, 2))(conv2)
    flat = layers.Flatten()(pool2)
    flat = layers.Dropout(0.5)(flat)
    fc = layers.Dense(256,activation = 'relu')(flat)
    fc2 = layers.Dense(2,activation = 'sigmoid')(fc)
    model = Model(inp,fc2)
    
    ##############################train and test##################
    model.compile(loss = 'categorical_crossentropy', optimizer = optimizers.Adam(lr = 1e-4), metrics = ['acc'])
    model.summary()   
    callbacks_list = [tf.keras.callbacks.ModelCheckpoint(filepath = 'earth.h5',monitor = 'val_loss',save_best_only = True), 
                      tf.keras.callbacks.ReduceLROnPlateau(monitor = 'val_loss', factor = 0.5, patience = 10),
                     ]#callbacks
    history = model.fit(x_train,y_train, validation_data=(x_test,y_test), epochs = 50, batch_size = 64,callbacks = callbacks_list) 
    
    ytrain_pred = model.predict([x_train])
    ytest_pred = model.predict([x_test])
    
    ytr_pred = ytrain_pred[:,0] > ytrain_pred[:,1]
    yte_pred = ytest_pred[:,0] > ytest_pred[:,1]
    
    tr = (ytr_pred == y_train[:,1])
    te = (yte_pred == y_test[:,1])
    hise.append(history.history)
    pree.append(precision_score(yte_pred,y_test[:,0]))
    rece.append(recall_score(yte_pred,y_test[:,0]))
    cone.append(confusion_matrix(yte_pred,y_test[:,0]))

res_hise = []
for i in range(10):
    res_hise.append(hise[i]['val_acc'][49])


fpe = 0
fne = 0
for i in range(len(cone)):
    fpe += cone[i][0,1]/len(cone)
    fne += cone[i][1,0]/len(cone)

a = 0
for i in range(10):
    a += np.array(hise[i]['val_acc'])/10
    
b = 0
for i in range(10):
    b += np.array(hisq[i]['val_acc'])/10  
    
c = 0
for i in range(10):
    c += np.array(hiso[i]['val_acc'])/10 
    
plt.plot(a)
plt.plot(b)
plt.plot(c)
    
        
    










