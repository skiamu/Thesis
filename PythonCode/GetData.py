#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jan 28 21:29:42 2018
script for retriving data from yahoo finance and get started with pandas
@author: andrea
"""
import numpy as np
import matplotlib.pyplot as plt
from sklearn import mixture
import scipy.io
# import data from the .mat file
Returns = scipy.io.loadmat('AssetClassReturns.mat') # it's a dictionary
Returns = Returns['Returns']

