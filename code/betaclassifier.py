#!/usr/bin/env python3

from math import log

from scipy.special import beta
import pandas as pd


class BetaClassifier:
    
    def __init__(self):
        
        pass
    
    
    def set_params(self, alpha0, alpha1, beta0, beta1, N0, N1):
        
        self.__alpha0 = alpha0    # alpha parameter of beta function for class 0
        self.__beta0 = beta0      # beta parameter of beta function for class 0
        self.__alpha1 = alpha1    # alpha parameter of beta function for class 1
        self.__beta1 = beta1      # beta parameter of beta function for class 1
        self.__N0 = N0            # number of subjects with negative class in training data
        self.__N1 = N1            # number of subjects with positive class in training data
        
        
        
    def predict_proba(self, X, uniq_CMV_label, uniq_total_label):
        
        ''' Predicts classes y for features X.
        Input:
        - X, pandas DataFrame
            X should be a pandas dataframe with features as columns and subjects as rows.
            The dataframe should contain two features:
                - # of unique CMV associated TCRs
                - total # of unique TCRs
        - uniq_CMV_label
            column name in X for # of unique CMV associated TCRs
        - uniq_total_label
            column name in X for total # of unique TCRs
        '''
        
        predictions = {'prediction_proba': []}
        
        for ix, row in X.iterrows():
            k = row[uniq_CMV_label]     # # of unique CMV associated TCRs
            n = row[uniq_total_label]   # total # of unique TCRs
            
            predictions['prediction_proba'].append(self._predict(k, n))
            
        return pd.DataFrame(predictions)
        
        
    def _predict(self, k, n):
        
        return log(self.__N1+1) - log(beta(self.__alpha1, self.__beta1)) + log(beta(k+self.__alpha1, n-k+self.__beta1)) - (
            log(self.__N0+1) - log(beta(self.__alpha0, self.__beta0)) + log(beta(k+self.__alpha0, n-k+self.__beta0))
        )