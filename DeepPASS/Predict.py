#!/usr/bin/env python
# coding: utf-8


'''

*Program:Predict.py
*Author:Guo-Hua yuan, Guo-Wei Li
*Function: predict sequence from or not from poly(A) site 
*Usage:
	python Predict.py -m model.h5 -p inputSeq.txt -o outputDir
		-m: DeepPASS model file
		-p: input sequence 
		(split by tab,
		first colum: seqID,
		second colum: sequence(200nt) )
		-o: outputdir
'''

import os
import argparse
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers, models
from tensorflow.compat.v1 import ConfigProto
from tensorflow.compat.v1 import InteractiveSession

#get filepath of input and output
def get_filepath():
    parser = argparse.ArgumentParser()
    parser.add_argument("-M","-m","--model", dest = "model", help = "Input model need to be evaluated")
    parser.add_argument("-P","-p","--predict",dest = "predict", help = "Set need to be predicted")
    parser.add_argument("-O","-o","--outputdir",dest = "outputdir",help = "Output dir of evaluation result")
    args = parser.parse_args()
    return args

# Convert data to tensor
def convert2tensor(data):
    NT_dict = {'A':np.array([1,0,0,0]),'T':np.array([0,1,0,0]),'U':np.array([0,1,0,0]),'G':np.array([0,0,1,0]),'C':np.array([0,0,0,1]),'N':np.array([0,0,0,0])}
    label_dict = {'0': 0, '1':1}
    
    seq = tf.convert_to_tensor(data["seq"].apply(lambda x:np.array([ NT_dict[i] for i in x],dtype=np.float32)))
    
    Set = {"seq":seq}
    
    return Set


if __name__ == "__main__":
    
    #################################### Configution ##################################
    # GPU Device
    if tf.test.is_gpu_available() == True:
        CONFIG = ConfigProto()
        CONFIG.gpu_options.allow_growth = True
        SESSION = InteractiveSession(config=CONFIG)


    ################################# Data Processing ################################
    print("\n********************** Start Data Processing **************************")
    # get input_file path
    args = get_filepath()
    model = args.model
    pre_f = args.predict
    output_dir = args.outputdir
    # load data
    dataset_pre = pd.read_csv(pre_f,sep='\t',index_col = False,header = None)

    # format:[seq]
    dataset = dataset_pre.iloc[:,-1:]
    dataset.columns = ["seq"]
    #convert to tensor
    PredictSet = convert2tensor(dataset)
    
    print("\n********************** Finished Data Processing ***********************")

    ################################### Prediction ####################################
    print("\n************************* Start Prediction ****************************")
    # Evaluate model performance on NewSet
    if not (os.path.exists(output_dir)):
        os.mkdir(output_dir)
    # Load best model
    os.environ["HDF5_USE_FILE_LOCKING"] = "FALSE"
    best_model =  tf.keras.models.load_model(os.path.join(model))
    
    # Predict
    prediction = best_model.predict(PredictSet["seq"])
    # Evaluate (negative: 0, positive: 1)
    pre_label = np.argmax(prediction,axis=1)
    
    # Output(original data with label)
    label_df = pd.DataFrame(pre_label,columns=["label"])
    total_df = pd.concat([dataset_pre,label_df],axis = 1)
    total_df.to_csv(os.path.join(output_dir,"Predict_Result.txt"),header = False,index = False,sep = '\t')
    
    
    print("\n********************** Finished Prediction ****************************")
    
    

