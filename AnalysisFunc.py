import numpy as np
from scipy.stats import ttest_ind, ttest_ind_from_stats
import matplotlib.pyplot as plt
from AdditionalFunc import*


# Data processing
def OTSU_enhance(img_gray, th_begin=0, th_end=256, th_step=1):
    assert img_gray.ndim <= 2, "must input a gary_img instead of {}d.".format(img_gray.ndim)
    
    if max(img_gray)!=0:
        img_gray = 255*img_gray/max(img_gray)
        
    max_g = 0
    suitable_th = 0
    for threshold in range(th_begin, th_end, th_step):
        bin_img = img_gray > threshold
        bin_img_inv = img_gray <= threshold
        fore_pix = np.sum(bin_img)
        back_pix = np.sum(bin_img_inv)
        if 0 == fore_pix:
            break
        if 0 == back_pix:
            continue
            
        w0 = float(fore_pix) / img_gray.size
        u0 = float(np.sum(img_gray * bin_img)) / fore_pix
        w1 = float(back_pix) / img_gray.size
        u1 = float(np.sum(img_gray * bin_img_inv)) / back_pix
        # intra-class variance
        g = w0 * w1 * (u0 - u1) * (u0 - u1)
        if g > max_g:
            max_g = g
            suitable_th = threshold
    a = img_gray[np.nonzero(img_gray>suitable_th)]
    b = img_gray[np.nonzero(img_gray<=suitable_th)]
    t, p = ttest_ind(a, b, equal_var=False)
    if p<0.05:
        diff = True
    else:
        diff = False
    diff = True
    return (suitable_th, diff, img_gray)

def PolarClassifier(profile, ratio=1/5, Mode="Polarity"):
    
    celllen = len(profile)
    threshold, diff, mod_profile = OTSU_enhance(profile)
    profile_2v = mod_profile>threshold # can be modified to real values.
    LPole = profile_2v[:int(celllen*ratio)]
    MCell = profile_2v[int(celllen*ratio):int(celllen*(1-ratio))]
    RPole = profile_2v[int(celllen*(1-ratio)):]
    
    L_verify = np.any(LPole)
    R_verify = np.any(RPole)
    M_verify = np.any(MCell)
    tmp = mod_profile.copy()
    mod_profile[np.nonzero(mod_profile<=threshold)]=0
    mod_profile[np.nonzero(mod_profile>threshold)]=1
    
#     fig = plt.figure(figsize=(8,3))
#     ax = fig.add_subplot(111)
#     ax.plot(np.arange(len(tmp)), tmp, color='g', linewidth=2)
    #ax.plot(np.arange(len(tmp)), np.ones(len(tmp))*threshold, color='grey', linewidth=2, linestyle='--')
#     ax.plot(np.arange(len(tmp)), np.ones(len(tmp))*np.mean(tmp), color='grey', linewidth=2, linestyle='--')
#     ax.plot(np.arange(len(tmp)), np.ones(len(tmp))*np.max(tmp), color='grey', linewidth=2, linestyle='--')
#     plt.grid('off')
#     plt.tight_layout()
#     ax.set_ylim([0, 255])
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#     ax.yaxis.set_ticks_position('left')
#     ax.xaxis.set_ticks_position('bottom')
#     ax.axhline(linewidth=2, y=0, color='k')
#     ax.axvline(linewidth=2, x=0, color='k')
#     ax.patch.set_facecolor('white')
#     plt.show()
    
    #PeakRange = (np.nonzero(mod_profile>0)[0][0], np.nonzero(mod_profile>0)[0][-1])
    
    if diff == False:
        PeakRange = np.arange(len(tmp))
        ind = PeakRange[tmp[PeakRange].argmax()]
        PolePos = 1 if ind/len(tmp)>0.5 else 0
        return ("diffuse", PeakRange, PolePos)
    elif L_verify==True and R_verify==True:
        if M_verify == True:
            PeakRange = np.arange(len(tmp))
            ind = PeakRange[tmp[PeakRange].argmax()]
            PolePos = 1 if ind/len(tmp)>0.5 else 0
            return ("diffuse", PeakRange, PolePos)
        else:
            PeakRange = np.nonzero(mod_profile>0)[0]
            ind = PeakRange[tmp[PeakRange].argmax()]
            PolePos = 1 if ind/len(tmp)>0.5 else 0
            return ("multipolar", PeakRange, PolePos)
    elif L_verify==False and R_verify==False:
        if M_verify == True:
            PeakRange = np.nonzero(mod_profile>0)[0]
            ind = PeakRange[tmp[PeakRange].argmax()]
            PolePos = 1 if ind/len(tmp)>0.5 else 0
            return ("midcell", PeakRange, PolePos)
        else:
            PeakRange = np.arange(len(tmp))
            ind = PeakRange[tmp[PeakRange].argmax()]
            PolePos = 1 if ind/len(tmp)>0.5 else 0
            return ("diffuse", PeakRange, PolePos)
    else:
        if M_verify == True:
            PeakRange = np.nonzero(mod_profile>0)[0]
            ind = PeakRange[tmp[PeakRange].argmax()]
            PolePos = 1 if ind/len(tmp)>0.5 else 0
            return ("occlusion", PeakRange, PolePos)
        else:
            PeakRange = np.nonzero(mod_profile>0)[0]
            ind = PeakRange[tmp[PeakRange].argmax()]
            PolePos = 1 if ind/len(tmp)>0.5 else 0
            return ("unipolar", PeakRange, PolePos)


def PolarityScore(data, ratio=1/3, PopZ=False):
    # "PopZ=0" means there is a PopZ focus at left side of a cell.
    # "PopZ=1" means there is a PopZ focus at right side of a cell.
    # "PopZ=3" means there is no PopZ focus in a cell.
    new_data = []
    length = []
    for i, d in enumerate(data):
        length.append(len(d))
        q1 = int(ratio*(len(d)))
        q3 = int((1-ratio)*(len(d)))
        PS = 0
        
        left = d[0:q1].max()
        right = d[q3:].max()
        
        if np.mean(d) == 0:
            PS = 0
        elif PopZ!=False:
            if PopZ[i]==0:
                PS = left/np.mean(d)
            else:
                PS = right/np.mean(d)
        else:
            if right > left:
                PS = right/np.mean(d)
            else:
                PS = left/np.mean(d)
        
        new_data.append(PS)
    
    return (np.array(length), np.array(new_data), np.mean(new_data), np.std(new_data))

        
def flip_profile(data, ratio=1/3, PopZ=False):
    new_data = []
    for i, d in enumerate(data):
        q1 = int(ratio*(len(d)))
        q3 = int((1-ratio)*(len(d)))
        try:
            left = d[0:q1].max()
            right = d[q3:].max()
        except:
            print('error here!!!', d)
            
        if type(PopZ)!=type(False):
            if PopZ[i]==0:
                new_data.append(d)
            else:
                reversed_arr = d[::-1]
                new_data.append(reversed_arr)
        else:
            if right > left:
                reversed_arr = d[::-1]
                new_data.append(reversed_arr)
            else:
                new_data.append(d)
    
    return np.array(new_data)


def arr_length_norm(arrs, method):
    min_len = 2000
    for arr in arrs:
        tmp_len = len(arr)
        if tmp_len < min_len:
            min_len = tmp_len
    
    if method == 'mean':
        transform_arrs = []
        for arr in arrs:
            tmp_len = len(arr)
            if tmp_len > min_len:
                tmp_arr = []
                size = tmp_len-min_len+1
                for i in range(min_len):
                    tmp_arr.append(np.mean(arr[i:i+size]))
                transform_arrs.append(np.array(tmp_arr))
            else:
                transform_arrs.append(arr.flatten())

        return np.array(transform_arrs)


def profile_statistics(data, Norm=True, NormV=None):
    EditProfile = []
        
    if Norm == True:
        for line in data:
            if np.max(line) == 0:
                line = np.array(line)/1
            else:
                line = np.array(line)/np.max(line)
            EditProfile.append(line)
    else:
        if NormV!=None:
            for line in data:
                line = np.array(line)/NormV
                EditProfile.append(line)
        else:
            EditProfile = data
        
    avg = np.mean(EditProfile, axis=0)
    std = np.std(EditProfile, axis=0)
    
    return (avg, std, EditProfile)


def TotalIntensity(data, mode="mean"):
    if mode == "mean":
        func = np.mean
    else:
        func = np.sum
    TI = []
    for line in data:
        TI.append(func(line))
    tot_avg = np.mean(TI)
    tot_std = np.std(TI)
    
    return (tot_avg, tot_std, TI)

        
def RegionCompare(profiles, ratio=1/3, PopZ=False): # need to be fixed (consider PopZ side)
    Pole, PoleMax = [], []
    Oppo, OppoMax = [], []
    Mid, MidMax = [], []
    length = []
    for i, profile in enumerate(profiles):
        # Record the length of the cell.
        length.append(len(profile))
        # Compute the index of the three regions.
        head = int(len(profile)*ratio)
        tail = int(len(profile)*(1-ratio))
        # Sum up the total intensity of the three regions.
        fist = np.sum(profile[0:head])
        second = np.sum(profile[head:tail])
        third = np.sum(profile[tail:])
        # Normalize the intensity by the total intensity.
        total = np.sum(fist)+np.sum(second)+np.sum(third)
        ft_norm = fist/total
        sd_norm = second/total
        td_norm = third/total
        # The maximum values in the three regions.
        max1 = np.max(profile[0:head])/np.mean(total)
        max2 = np.max(profile[head:tail])/np.mean(total)
        max3 = np.max(profile[tail:])/np.mean(total)
        # Flip the intensity which locates at the PopZ to the left side of the cell.
        if type(PopZ)!=type(False):
            if PopZ[i]==0:
                Pole.append(ft_norm)
                Mid.append(sd_norm)
                Oppo.append(td_norm)
                PoleMax.append(max1)
                MidMax.append(max2)
                OppoMax.append(max3)
            else:
                Pole.append(td_norm)
                Mid.append(sd_norm)
                Oppo.append(ft_norm)
                PoleMax.append(max3)
                MidMax.append(max2)
                OppoMax.append(max1)
        else:
            if ft_norm > td_norm:
                PoleMax.append(max1)
                MidMax.append(max2)
                OppoMax.append(max3)
                Pole.append(ft_norm)
                Mid.append(sd_norm)
                Oppo.append(td_norm)
            else:
                Pole.append(td_norm)
                Mid.append(sd_norm)
                Oppo.append(ft_norm)
                PoleMax.append(max3)
                MidMax.append(max2)
                OppoMax.append(max1)
    # Statistics
    Pole_avg, Mid_avg, Oppo_avg = np.mean(Pole), np.mean(Mid), np.mean(Oppo)
    Pole_std, Mid_std, Oppo_std = np.std(Pole), np.std(Mid), np.std(Oppo)
    PoleM_avg, MidM_avg, OppoM_avg = np.mean(PoleMax), np.mean(MidMax), np.mean(OppoMax)
    PoleM_std, MidM_std, OppoM_std = np.std(PoleMax), np.std(MidMax), np.std(OppoMax)
    # Normalize the composition ratio.
    scaler = 1/(Pole_avg+Mid_avg+Oppo_avg)
    avg = (scaler*Pole_avg, scaler*Mid_avg, scaler*Oppo_avg)
    std = (scaler*Pole_std, scaler*Mid_std, scaler*Oppo_std)
    avgM = (PoleM_avg, MidM_avg, OppoM_avg)
    stdM = (PoleM_std, MidM_std, OppoM_std)
    return (avg, std, avgM, stdM)
    

def MidIntensity(profile, ratio):
    midregion_head = int(len(profile)*ratio)
    midregion_tail = int(len(profile)*(1-ratio))
    midtotal = np.sum(profile[midregion_head:midregion_tail])
    total = np.sum(profile)
    return midtotal/total


def LenNorm(profile, Len=150):
    profile = profile.flatten()
            
    def linearEq(x1, y1, x2, y2, value):
        a = (y2-y1)/(x2-x1)
        b = y1-a*x1
        return a*value + b
    
    def area(x1, y1, x2, y2):
        return (y1+y2)*abs(x2-x1)/2
    
    
    y_collect = []
    Grid = Len-1
    arr = np.zeros(Grid)
    profile_ind = np.arange(len(profile))
    if len(profile) != Len:
        stride = (len(profile)-1)/(Len-1)
        for i in range(Grid):
            if i == Grid-1:
                x_last = i*stride
                y_last = y_this
                ceil = np.nonzero(x_last<=profile_ind)[0][0]
                x_this = profile_ind[-1]
                y_this = profile[x_this]
                midpoints = x_this - ceil
                y_collect.append(y_this)
                
                if midpoints == 0:
                    arr[i] = area(x_last, y_last, x_this, y_this)
                else:
                    cover = 0
                    for j in range(midpoints+1):
                            if j==0:
                                cover+=area(x_last, y_last, ceil+j, profile[ceil+j])
                            elif j==midpoints:
                                cover+=area(ceil+j-1, profile[ceil+j-1], x_this, y_this)
                            else:
                                cover+=area(ceil+j-1, profile[ceil+j-1], ceil+j, profile[ceil+j])
                    arr[i] = cover
            else:
                x_this = (i+1)*stride
                x_last = i*stride
                ceil = np.nonzero(x_this<=profile_ind)[0][0]
                floor = np.nonzero(x_last<=profile_ind)[0][0]
                midpoints = ceil-floor
                x2, x1 = ceil, ceil-1
                y2, y1 = profile[x2], profile[x1]
                y_this = linearEq(x1, y1, x2, y2, x_this)


                if i==0:
                    x0, y0 = 0, profile[0]
                    y_collect.append(y0)
                    if midpoints == 0:
                        arr[0] = area(x0, y0, x_this, y_this)
                    else:
                        cover = 0
                        for j in range(midpoints+1):
                            if j==0:
                                cover+=area(x0, y0, floor+j, profile[floor+j])
                            elif j==midpoints:
                                cover+=area(floor+j-1, profile[floor+j-1], x_this, y_this)
                            else:
                                cover+=area(floor+j-1, profile[floor+j-1], floor+j, profile[floor+j])
                        arr[0] = cover
                else:
                    if midpoints == 0:
                        arr[i] = area(x_last, y_last, x_this, y_this)
                    else:
                        cover = 0
                        for j in range(midpoints+1):
                            if j==0:
                                cover+=area(x_last, y_last, floor+j, profile[floor+j])
                            elif j==midpoints:
                                cover+=area(floor+j-1, profile[floor+j-1], x_this, y_this)
                            else:
                                cover+=area(floor+j-1, profile[floor+j-1], floor+j, profile[floor+j])
                        arr[i] = cover

                y_last = y_this
                y_collect.append(y_this)

    else:
        y_collect.append(profile[0])
        for i in range(Len-1):
            x2, x1 = profile_ind[i+1], profile_ind[i]
            y2, y1 = profile[x2], profile[x1]
            y_collect.append(y2)
            cover = area(x1, y1, x2, y2)
            arr[i] = cover

    return arr, y_collect

def Batch_LenNorm(profiles, option='points', Length=150):
    res_collect = []
    for profile in profiles:
        arr, y_collect = LenNorm(profile, Len=150)
        if option=="area":
            res_collect.append(arr)
        else:
            res_collect.append(y_collect)
            
    return res_collect


def MassCenterScore(data, PopZ=False):
    new_data = []
    length = []
    for i, line in enumerate(data):
        length.append(len(line))
        if MassCenter1D(line)>len(line)-1:
            print('bugbug')
        a = 1-(MassCenter1D(line)/(len(line)-1))
        b = (MassCenter1D(line)/(len(line)-1))
        if type(PopZ)!=type(False):
            if PopZ[i]==0:
                new_data.append(b)
            else:
                new_data.append(a)
        else:
            if a>b:
                new_data.append(b)
            else:
                new_data.append(a)
        
    return (np.array(length), np.array(new_data), np.mean(new_data), np.std(new_data))


def PolarLocalRatio(data, ratio=1/3, PopZ=False, PeakRange=[]):
    new_data = []
    length = []
    for i, d in enumerate(data):
        length.append(len(d))
        q1 = int(ratio*(len(d)))
        q3 = int((1-ratio)*(len(d)))
        tot = np.sum(d)
        PUL = 0
        
        left = np.sum(d[0:q1])
        right = np.sum(d[q3:])
        if tot == 0:
            PUL = 0
        else:
            if type(PopZ)!=type(False):
                if PeakRange:
                    PUL = np.sum(d[PeakRange[i]])/tot
                else:
                    if PopZ[i]==0:
                        PUL = left/tot
                    else:
                        PUL = right/tot
            else:
                if right > left:
                    PUL = right/tot
                else:
                    PUL = left/tot
        
        new_data.append(PUL)
    
    return (np.array(length), np.array(new_data), np.mean(new_data), np.std(new_data))

def PatternNormalize(profiles1, profiles2):
    # profiles1 is the reference of the average spatial pattern used to normalize.
    # profiles1 and 2 have to have the same length (dimension).
    profile_collect = []
    for profile in profiles2:
        profile_collect.append(profile/profiles1)
    return profile_collect


def Polar2ndChannel(data2, PeakRange, mode='colocalization'):
    """
    !!!The profiles in the data should have normalized length and correct direction.
    """
    new_data = []
    for i, d in enumerate(zip(data2, PeakRange)):
        ind_arr = np.array(d[1])
        sig_arr = np.array(d[0])
        polar_tot = np.sum(sig_arr[ind_arr])
        diff_pixel = (np.sum(sig_arr)-polar_tot)/(len(sig_arr)-len(ind_arr))
        B = len(sig_arr)*diff_pixel
        A = polar_tot-(diff_pixel*len(ind_arr))
        if A+B == 0:
            Bowman_colocal = 0
        else:
            Bowman_colocal = A/(A+B)
        if mode == 'colocalization':
            new_data.append(Bowman_colocal)
        else:
            new_data.append(polar_tot)
        
    return (np.array(new_data), np.mean(new_data), np.std(new_data))


def LR_diff(data, PopZ=False, ratio=1/3, pseudo=0):
    new_data = []
    length = []
    for i, line in enumerate(data):
        q1 = int(ratio*(len(line)))
        q3 = int((1-ratio)*(len(line)))
        tmpL = len(line[:q1])
        L = np.sum(line[:q1])+pseudo
        R = np.sum(line[-tmpL:])+pseudo
        if len(line[:q1])!=len(line[-tmpL:]):
            print('len error')
            exit()
        if L==0 or R==0:
            #new_data.append(1)
            print('skip')
        else:
            if PopZ!=False:
            #if type(PopZ)!=type(False):
                if PopZ[i]==0:
                    new_data.append(L/R)
                else:
                    new_data.append(R/L)
            else:
                if (L/R)>(R/L):
                    new_data.append(L/R)
                else:
                    new_data.append(R/L)    
    return (np.array(length), np.array(new_data), np.mean(new_data), np.std(new_data))


def LR_max(data, PopZ=False, ratio=1/3, pseudo=10**(-5)):
    new_data = []
    length = []
    for i, line in enumerate(data):
        q1 = int(ratio*(len(line)))
        q3 = int((1-ratio)*(len(line)))
        L = np.max(line[:q1])+pseudo
        R = np.max(line[q3:])+pseudo

        a = L/R
        b = R/L
        if PopZ!=False:
        #if type(PopZ)!=type(False):
            if PopZ[i]==0:
                new_data.append(a)
            else:
                new_data.append(b)
        else:
            if a>b:
                new_data.append(a)
            else:
                new_data.append(b)
        
    return (np.array(length), np.array(new_data), np.mean(new_data), np.std(new_data))

def Pearson_Correlation(x_arrs, y_arrs, ratio=1/3, pseudo_v=10**(-5)):
    # pearson coeff = cov(x, y) / ( sig(x) * sig(y) )
    res = []
    for x_arr, y_arr in zip(x_arrs, y_arrs):
        r = np.corrcoef(x_arr.flatten(), y_arr.flatten())
        if np.isnan(r[0][1]):
            pass
#             print('zero std bug', r, np.std(x_arr), np.std(y_arr))
        else:
            res.append(r[0][1])
    
    return res


def PeakFeature(profile, PeakRange):
    height = np.max(profile)
    width = len(PeakRange)
    if height<0:
        print('width bug')
    return height/width


def kymo_preprocess(profileCh1, profileCh2):
    """
    ::profileCh1::
    The input is an one-layer array.
    It contains [[mother cell array], [daughter cell array]].
    Besides, the second input ::profileCh2:: represents the second channel signals.
    """
    # Always put the parental cell (with higher PopZ peak) at left side.
#     if np.any(profileCh1[1]): # If and only if the second cell does not contain nothing.
#         if max(profileCh1[0])<max(profileCh1[1]):
#             # Exchange the first channel of two cell arrays.
#             tmp = profileCh1[0]
#             profileCh1[0] = profileCh1[1]
#             profileCh1[1] = tmp
#             # Exchange the second channel of the two cell arrays.
#             tmp = profileCh2[0]
#             profileCh2[0] = profileCh2[1]
#             profileCh2[1] = tmp

    # Always put the PopZ peak on left-hand side.
    if profileCh1[0].argmax()/len(profileCh1[0]) >= 0.5:
        profileCh1[0] = profileCh1[0][::-1]
        profileCh2[0] = profileCh2[0][::-1]

    # Merge the parental and daughter cell profiles.
    if np.any(profileCh1[1]): # If and only if the second cell does not contain nothing.
        kymo_tmp1 = np.concatenate((profileCh1[0].flatten(), profileCh1[1].flatten()))
        kymo_tmp2 = np.concatenate((profileCh2[0].flatten(), profileCh2[1].flatten()))
        if np.max(profileCh2[1].flatten()) == 0:
            daughter = 0
        else:
            #sum
            daughter = np.sum(profileCh2[1].flatten())#np.max(profileCh2[1].flatten()))
        if np.max(profileCh2[0].flatten()) == 0:
            mother = 0
        else:
            #sum
            mother = np.sum(profileCh2[0].flatten())#np.max(profileCh2[0].flatten()))
        
    else: # There is only one cell.
        kymo_tmp1 = profileCh1[0].flatten()
        kymo_tmp2 = profileCh2[0].flatten()
        if np.max(profileCh2[0].flatten()) == 0:
            mother = 0
        else:
            #sum
            mother = np.sum(profileCh2[0].flatten())
        daughter = 0
    kymo1 = (kymo_tmp1)#/kymo_tmp1.max())
    kymo2 = (kymo_tmp2)#/kymo_tmp2.max())

    if np.isnan(mother) or np.isnan(daughter):
        print('nan error')
    return (kymo1, kymo2, mother, daughter)