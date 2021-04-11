from scipy.io import loadmat
import numpy as np
import csv
import os
import sys
from collections import OrderedDict


def read_mat(parent_path, filename, middle_path=None):
    pathProg = parent_path
    os.chdir(pathProg)
    package_type = "loadmat"
    if type(middle_path)==type('abc'):
        os.chdir(parent_path)
        tmp = parent_path + r'/{}'.format(middle_path)
        os.chdir(tmp)
    try:
        # load .mat file by scipy.io.loadmat (usually for old .mat)
#         print('Prepare to process {}'.format(filename))
        f = loadmat(filename) # e.g. ("SPOT_OBJ_RFP.mat")
        data = f['cellList']['meshData']
        package_type = "loadmat"
    except:
        # load .mat file by h5py (usually for new .mat)
        print("The file is broken or not completed.")
        data = None
        exit()
    return (package_type, data)

def preprocess(file, Cells, signal, TwoInOne=False):
    ID_start = Cells['count'] # The last ID in dictionary.
    cells_file = file[0][0][0][0][0] # the arribute of the signal
    cell_num = len(cells_file) # the number of cells
    #print('ID_start', ID_start)
    print('Processing '+str(cell_num)+' cells in signal'+str(signal)+'.')
    adjustment = 0
    for i in range(cell_num):
        # Read info of i-th cell
        Base = cells_file[i][0][0]
        if len(Base['signal1'])<3 or len(Base['signal1'])>150:
            print("Warning: a cell with length {} shorter than 3 or longer than 150 units will be skipped.".format(len(Base['signal1'])))
            adjustment += 1
            continue
        # Adjust ID number
        ID = ID_start+i+1-adjustment
        # Create a dictionary for recording "this" cell.
        Cells['{}'.format(ID)] = Cells['{}'.format(ID)] if Cells.get('{}'.format(ID), 0) else {}
        thisCell = Cells['{}'.format(ID)]
        # Create a subdictionary for recording "the signal".
        thisCell['S{}'.format(signal)] = {}
        thisSignal = thisCell['S{}'.format(signal)]
        if TwoInOne == True:
            # Save FP profile.
            thisSignal['profile'] = Base['signal{}'.format(signal+1)]
        else:
            # Save FP profile.
            thisSignal['profile'] = Base['signal1']
        update_ID = ID
    return update_ID, Cells

def grouping(parent_path, exp, labels, ChNum=1):
    """
    ::parent_path:: 
    [Original type of string] The path for finding the files.
    ::exp:: 
    [Numpy array] A deep structure containing all groups of data.
    ::labels:: 
    [Numpy array] Names for all groups of data.
    ::TwoInOne::
    [Integer] Indicate how many FP channels included in one file.
    """
    res = OrderedDict()
    for label, group in zip(labels, exp): # group1, group2...
        print("Start to process the group of", label, ".")
        res['{}'.format(label)] = {}
        Cells = res['{}'.format(label)]
        Cells['count'] = -1
        Cells['PopZ'] = 1
        switch = 0
        for img in group: # img1, img2...
            if "TwoInOne" in img:
                switch = 1
            elif "NP" in img:
                Cells['PopZ'] = 0
            else:
                if switch:
                    for signal in range(ChNum): # e.g. file = RFP file
                        # Open the file
#                         print("Open and save the file", img[0], "to the temporary space.")
                        check, data = read_mat(parent_path, img[0])
                        # Save the data in this file into a dictionary
                        count_update, Cells = preprocess(data, Cells, signal, TwoInOne=True)
                else:
                    for signal, file in enumerate(img): # e.g. file = RFP file
                        # Open the file
#                         print("Open and save the file", file, "to the temporary space.")
                        check, data = read_mat(parent_path, file)
                        # Save the data in this file into a dictionary
                        count_update, Cells = preprocess(data, Cells, signal)

                Cells['count'] = count_update
            
    
    return res


def TL_preprocess(file, Cells, signal):
    ID = Cells['count'] # The last ID in dictionary.
    cells_file = file[0][0][0] # the arribute of the signal
    cell_num = len(cells_file) # the number of cells
#     print('Processing '+str(cell_num)+' cells in signal'+str(signal)+'.')
    # Create a dictionary for recording "this" cell.
    Cells['{}'.format(ID)] = Cells['{}'.format(ID)] if Cells.get('{}'.format(ID), 0) else {}
    thisCell = Cells['{}'.format(ID)]
    for i in range(cell_num):
        # Read info of i-th cell
        time_base = cells_file[i][0]
        for j in range(time_base.shape[0]):
            Base = time_base[j][0][0]

            if len(Base['signal1'])<3:
                print("Warning: a cell shorter than 3 units will be skipped.")
                continue
            # Adjust time series.
            # Create a subdictionary for recording "the signal".
            thisCell['T{}'.format(i)] = thisCell['T{}'.format(i)] if thisCell.get('T{}'.format(i), 0) else {}
            thisCell['T{}'.format(i)]['L{}'.format(j)] = thisCell['T{}'.format(i)]['L{}'.format(j)] if thisCell['T{}'.format(i)].get('L{}'.format(j), 0) else {}
            thisCell['T{}'.format(i)]['L{}'.format(j)]['S{}'.format(signal)] = {}
            thisSignal = thisCell['T{}'.format(i)]['L{}'.format(j)]['S{}'.format(signal)]
            # Save FP profile.
            thisSignal['profile'] = Base['signal{}'.format(signal+1)]

    update_ID = ID+1
    return update_ID, Cells



def TL_preprocess2(file, Cells, signal, ID, t):
    cells_file = file[0][0][0][0][0] # the arribute of the signal
#     print('Processing cells in signal'+str(signal)+'.')
    thisCell = Cells['{}'.format(ID)]
    for j in range(len(cells_file)):
        Base = cells_file[j]
        # Create a subdictionary for recording "the signal".
        thisCell['T{}'.format(t)] = thisCell['T{}'.format(t)] if thisCell.get('T{}'.format(t), 0) else {}
        thisCell['T{}'.format(t)]['L{}'.format(j)] = thisCell['T{}'.format(t)]['L{}'.format(j)] if thisCell['T{}'.format(t)].get('L{}'.format(j), 0) else {}
        thisCell['T{}'.format(t)]['L{}'.format(j)]['S{}'.format(signal)] = {}
        thisSignal = thisCell['T{}'.format(t)]['L{}'.format(j)]['S{}'.format(signal)]
        # Save FP profile.
        thisSignal['profile'] = Base['signal{}'.format(signal+1)]

    return Cells


def SingleCellTimeLapse(parent_path, exp, TwoInOne=2):
    """
    ***Input -exp- Structure***
    --exp: includ experimental and control group.
    ----group: contain all img (every one cell).
    ------img: contain all channel files (here should be one file.)
    --------file: .mat file generated from Oufti.
    
    ::parent_path:: 
    [Original type of string] The path for finding the files.
    ::exp:: 
    [Numpy array] A deep structure containing all groups of data.
    ::labels:: 
    [Numpy array] Names for all groups of data.
    ::TwoInOne::
    [Integer] Indicate how many FP channels included in one file.
    """
    res = OrderedDict()
    for group_key in exp: # group1, group2...
#         print("Start to process the group of", group_key, ".")
        res['{}'.format(group_key)] = {}
        Cells = res['{}'.format(group_key)]
        Cells['count'] = 0
        Cells['PopZ'] = 1
        for ID in exp['{}'.format(group_key)].keys(): # 1, 2...
            Cells['{}'.format(ID)] = Cells['{}'.format(ID)] if Cells.get('{}'.format(ID), 0) else {}
            for t, img in enumerate(exp['{}'.format(group_key)]['{}'.format(ID)]):# img1, img2...
                if "NP" in img:
                    Cells['PopZ'] = 0
                else:
                    if TwoInOne > 0:
                        for signal in range(TwoInOne): # e.g. file = RFP file
                            # Open the file
                            #print("Open and save the file", img, "to the temporary space.")
                            for pos, char in enumerate(img):
                                #if char == 's':
                                if char=='-':
                                    start = pos
                            middle_path = img[:start]
                            middle_path = '('+middle_path+')'
                            #print('middle_path', middle_path)
                            check, data = read_mat(parent_path, img, middle_path=middle_path)
                            # Save the data in this file into a dictionary
                            Cells = TL_preprocess2(data, Cells, signal, ID, t)
                    else:
                        print('The number of signal should be more than 0.')
    
    return res