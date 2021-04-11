import numpy as np
from collections import OrderedDict
import os, errno, shutil, re

def InputFromTXT(txt_path=''):
    if len(txt_path)!=0:
        path = txt_path
    else:
        path = r'./DataInput.txt'
    
    Data = {}
    ctrl = 0
    mode = ''
    with open(path, 'r') as f:
        for row in f:
            if "/*PLEASE DO NOT EDIT THE ABOVE REGION*/" in row:
                ctrl=1
            if ctrl==1:
                if 'parent_path' in row:
                    Data['parent_path'] = row.split(':')[1].strip()
                elif 'input_mode' in row:
                    mode = row.split(':')[1].strip().lower()
                elif 'group_number' in row:
                    Data['exp'] = [[] for _ in range(int(row.split(':')[1]))]
                elif 'channels' in row:
                    values = row.split(':')[1].strip()
                    value = values.split(';')
                    for ele, v in zip(Data['exp'], value):
                        if int(v)==2:
                            ele.append(np.array(['TwoInOne']))
                elif 'PopZ' in row:
                    values = row.split(':')[1]
                    value = values.split(';')
                    for ele, v in zip(Data['exp'], value):
                        if v.upper().strip()=='N':
                            ele.append(np.array(['NP']))
                elif 'labels' in row:
                    values = row.split(':')[1].strip()
                    value = values.split(';')
                    tmp = []
                    for v in value:
                        tmp.append(v)
                    Data['labels'] = np.array(tmp)
                elif 'name' in row:
                    label = row.split(':')[0]
                    values = row.split(':')[1].strip()
                    
                    if 'file' in mode.lower():
                        value = values.split(';')
                        tmp_f = []
                        for v in value:
                            tmp_f.append(v)
                        Data['exp'][int(re.findall('\d+', label)[0])-1].append(np.array(tmp_f))
                    elif 'folder' in mode.lower():
                        f_path = Data['parent_path'] + r'/{}'.format(values)
                        m_files = os.listdir(f_path)
                        tmp_f = [ff for ff in m_files if '.mat' in ff]
                        Data['exp'][int(re.findall('\d+', label)[0])-1].append(np.array(tmp_f))
                    
    Data['exp'] = np.array(Data['exp'])
    
    return Data



"""
Data for plots in paper.
"""
def fig1():
    parent_path = r'/home/r04b43015/Download/TEV/fig1/'
    group1, group2, group3 = [], [], []
    
    for mark in ['02', '03', '04']:
        img = np.array(['116DivIVA-GFP-{0}-Image-Export-{1}.mat'.format(mark, mark)])
        group1.append(img)
    
    for mark in ['02', '03']:
        img = np.array(['116DivIVA-GFP-M0052-{0}-Image-Export-{1}.mat'.format(mark, mark)])
        group2.append(img)
        
    for mark in ['02', '03', '05']:
        img = np.array(['116DivIVA-GFP-tevR-M0050-{0}-Image-Export-{1}.mat'.format(mark, mark)])
        group3.append(img)
    
    
    
    exp = np.array([group1, group2, group3])
    labels = np.array(['116DG', '116DG52', '116DG50'])
    Data = {}
    Data['parent_path'] = parent_path
    Data['exp'] = exp
    Data['labels'] = labels
    
    return Data


def fig2():
    parent_path = r'/home/r04b43015/Download/TEV/fig2'
    group1, group2, group3 = [], [], []
    
    for mark in ['02', '03']:
        img = np.array(['116DivIVA-GFP-M0050-{}.mat'.format(mark)])
        group1.append(img)
    
    for mark in ['08', '07']:
        img = np.array(['106TEV-116DivIVA-GFP-tevR-M0050-{}.mat'.format(mark)])
        group3.append(img)
    
    for mark in ['02', '03', '05']:
        img = np.array(['116DivIVA-GFP-tevR-M0050-{}.mat'.format(mark)])
        group2.append(img)
    
    
    exp = np.array([group1, group2, group3])
    labels = np.array(['DG50', 'DGS50', 'TDGS50'])
    Data = {}
    Data['parent_path'] = parent_path
    Data['exp'] = exp
    Data['labels'] = labels
    
    return Data


def fig3():
    parent_path = r'/home/r04b43015/Download/TEV/fig3'
    group1, group2, group3 = [np.array(['TwoInOne'])],[np.array(['TwoInOne'])],[np.array(['TwoInOne'])]
    
    for mark in ['03', '04', '05', '06', '07', '08']:
        img = np.array(['180302-1-{}.mat'.format(mark)])
        group1.append(img)
    
    for mark in ['02', '03', '04', '05', '06', '08', '07']:
        img = np.array(['180207-5-{}.mat'.format(mark)])
        group2.append(img)
    
    for mark in ['02', '03', '04']:
        img = np.array(['180207-4-{}.mat'.format(mark)])
        group3.append(img)
    
    
    exp = np.array([group1, group2, group3])
    labels = np.array(['PS', 'PSNC', 'PSNSSC'])
    Data = {}
    Data['parent_path'] = parent_path
    Data['exp'] = exp
    Data['labels'] = labels
    
    return Data

def fig4():
    parent_path = r'/home/r04b43015/Download/TEV/fig4'
    group1, group2 = [np.array(['TwoInOne'])], [np.array(['TwoInOne']),np.array(['NP'])]
    
    for mark in ['07', '03', '04']:
        img = np.array(['116RP-TEV(NSSC)-116DG(tevR-M0050)(50uM.1hr)-{}.mat'.format(mark)])
        group1.append(img)
    
    for mark in ['01', '04', '03', '05']:
        img = np.array(['TEV(NSSC)-116DG(tevR-M0050)(50uM.1hr)-{}.mat'.format(mark)])
        group2.append(img)
    
    
    exp = np.array([group1, group2])
    labels = np.array(['NSSCP', 'NSSC'])
    Data = {}
    Data['parent_path'] = parent_path
    Data['exp'] = exp
    Data['labels'] = labels
    
    return Data



def fig2_37c():
    parent_path = r'/home/r04b43015/Download/TEV/fig2(37c)'
    group1, group2, group3 = [], [], []
    
    for mark in ['04']:
        img = np.array(['116DG(m0050)(37c)-02-Image-Export-{}.mat'.format(mark)])
        group1.append(img)
    
    for mark in ['02']:
        img = np.array(['106TEV-116DG(tevr-m0050)(37c)-02-Image-Export-{}.mat'.format(mark)])
        group3.append(img)
    
    for mark in ['06', '07']:
        img = np.array(['116DG(tevR-m0050)(37c)-Image-Export-{}.mat'.format(mark)])
        group2.append(img)
    
    
    exp = np.array([group1, group2, group3])
    labels = np.array(['DG50', 'DGS50', 'TDGS50'])
    Data = {}
    Data['parent_path'] = parent_path
    Data['exp'] = exp
    Data['labels'] = labels
    
    return Data