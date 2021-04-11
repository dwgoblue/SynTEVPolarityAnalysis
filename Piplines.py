import numpy as np
from AnalysisFunc import*
from ReadDataFunc import grouping, SingleCellTimeLapse
from PlotFunc import PlotMethods
from DataInput import*
import pandas as pd
import warnings

class PiplinesObj():
    
    warnings.filterwarnings("ignore")
    
    def __init__(self, DataSet, columns=[], TimeLapse=False,
                 ChNum=2, mode='func', txt_path=r'{}/DataInput.txt'.format(os.getcwd())):
        if mode=='txt':
            r = InputFromTXT(txt_path=txt_path)
        else:
            r = DataSet()
        if TimeLapse==False:
            Data = grouping(r['parent_path'], r['exp'], r['labels'], ChNum)
                  
        else:
            Data = SingleCellTimeLapse(r['parent_path'], r['exp'], ChNum)
        self.ChNum = ChNum
        self.Data = Data
        self.PltFunc = PlotMethods(r['parent_path'], columns)

    def PopZOnlyPipline(self, plots):
        data = self.Data
        # Create a dictionary for storage.
        ProcessedData = OrderedDict()
        # Create temporary space for data storage.
        polarity_data = []
        PS_collect = []
        name_collect = []
        unipolar_collect = []
        feature_collect = []
        Len_collect = []
        TotInt_Collect = []
        LR_collect = []
        # Create a iteration objects that are indexing by the plasmids.
        item = list(data.items())
        for subfolder in item:
            # Initialize the space for saving all info.
            polarity_collect = []
            unipolar = 0
            diffuse = 0
            others = 0
            unipolar_group = np.zeros(3)
            profile_collect = []
            peak_feature = []
            Len_tmp = []
            CH1_index = []
            # The length of cell list.
            IDlen = (len(subfolder[1])-2)/3
            print('The group ', subfolder[0], 'has n=', IDlen*3)
            for i in range(len(subfolder[1])-2):
                profile = subfolder[1]['{}'.format(i)]['S0']['profile']
                profile_collect.append(profile)
                Len_tmp.append(len(profile))
                # Classify the spatial pattern of cells.
                polarity = PolarClassifier(profile, 1/3)
                if i/len(subfolder[1])-2 > 1/3:
                    gp_ctrl = 1
                elif i/len(subfolder[1])-2 > 2/3:
                    gp_ctrl = 2
                else:
                    gp_ctrl = 0
                    
                if polarity[0] == 'unipolar':
                    unipolar += 1
                    unipolar_group[gp_ctrl] += 1
                elif polarity[0] == 'diffuse':
                    diffuse += 1
                else:
                    others += 1
                peak_feature.append(PeakFeature(profile, polarity[1]))
                CH1_index.append(0) # 0: PopZ pole at left
            if len(subfolder[1])-2>180:
                sampnum = np.random.choice(np.arange(130, 170))
                select_ind = np.random.choice(np.arange(len(subfolder[1])-2), sampnum, replace=False)
            else:
                select_ind = np.arange(len(subfolder[1])-2)
            select_ind = select_ind.astype(int)
            with open('SampleNumber.txt', 'a+') as f:
                f.write('{0} n={1} \n'.format(subfolder[0], len(select_ind)))
                
            tot = TotalIntensity([profile_collect[ind] for ind in select_ind], mode="sum")
            TotInt_Collect.append(tot)
            Len_collect.append(Len_tmp)
            feature_collect.append([peak_feature[ind] for ind in select_ind])
            # Collect the name of the plasmid that transformed the cell.
            name_collect.append(subfolder[0])
            # Match the plasmid name and the cell number of its corresponding polarity
            unipolar_collect.append((np.mean(unipolar_group/IDlen),
                                     np.std(unipolar_group/IDlen),
                                     (unipolar_group/IDlen)))
            polarity_data.append((subfolder[0], (unipolar, diffuse, others)))
            # Draw profiles of PopZ FP distribution.
            #--> Calculate Polarity Score and cell length
            PS = PolarityScore([profile_collect[ind] for ind in select_ind], PopZ=CH1_index)
            PS_collect.append((PS[0], PS[1]))
            #--> Process the profiles
            EditProfiles = flip_profile([profile_collect[ind] for ind in select_ind], 1/3, PopZ=CH1_index)
            EditProfiles = Batch_LenNorm(EditProfiles, option='points')
            avg_arr, std_arr, EditProfiles = profile_statistics(EditProfiles, Norm=True)
            #--> Single cell traces (marked one line).
            self.PltFunc.ProfilePlot(EditProfiles, subfolder[0], mean_line=avg_arr)
            #--> Average of profiles with STD.
            self.PltFunc.ProfileStatPlot(avg_arr, std_arr, subfolder[0])
            LR_collect.append(LR_diff(profile_collect, ratio=1/3, PopZ=CH1_index))
            
        # Draw ratio of different kinds of polarity.
        max_bound = []
        for arr in TotInt_Collect:
            arr = arr[2]
            Q1 = np.percentile(arr, 75)
            Q3 = np.percentile(arr, 25)
            max_bound.append(Q1+1.5*(Q1-Q3))
        max_bound = max(max_bound)
        
        if 'tot' in plots:
            self.PltFunc.ViolinPlot([v[2]/max_bound for v in TotInt_Collect],
                                 name_collect,
                                 'TotalIntensity')
            df = pd.DataFrame(data=[v[2]/max_bound for v in TotInt_Collect], index=name_collect)
            df.to_csv('TotalIntensity.csv')
            
        if 'lr' in plots:
            self.PltFunc.ViolinPlot([v[1] for v in LR_collect],
                                 name_collect,
                                'LRratio')
            df = pd.DataFrame(data=[v[1] for v in LR_collect], index=name_collect)
            df.to_csv('TotalIntensity.csv')
        
        # Draw polarity analysis
        #--> Polarity Score with bar plot
        #--> Peak Score
        if 'ps' in plots:
            self.PltFunc.ViolinPlot([v[1] for v in PS_collect],
                                 name_collect,
                                 'PolarityScore')
            df = pd.DataFrame(data=[v[1] for v in PS_collect], index=name_collect)
            df.to_csv('PolarityScore.csv')
        
        if 'pks' in plots:
            self.PltFunc.ViolinPlot([np.array(arr) for arr in feature_collect],
                                 name_collect,
                                 'PeakScore')
            df = pd.DataFrame(data=[v for v in feature_collect], index=name_collect)
            df.to_csv('PeakScore.csv')
        
        print('The analysis is finishing.')


    def DiffuseSignalPipline(self, plots):
        data = self.Data
        # Draw variations of total FP intensity.
        name_collect = []
        TotInt_Collect = []
        item = list(data.items())
        for subfolder in item:
            profile_collect = []
            for i in range(len(subfolder[1])-2):
                profile = subfolder[1]['{}'.format(i)]['S0']['profile']
                profile_collect.append(profile)
            if len(subfolder[1])-2>180:
                sampnum = np.random.choice(np.arange(130, 170))
                select_ind = np.random.choice(np.arange(len(subfolder[1])-2), sampnum, replace=False)
            else:
                select_ind = np.arange(len(subfolder[1])-2)
            select_ind = select_ind.astype(int)
            with open('SampleNumber.txt', 'a+') as f:
                f.write('{0} n={1} \n'.format(subfolder[0], len(select_ind)))
            name_collect.append(subfolder[0])
            # Compute the statistics of total FP intensity.
            tot = TotalIntensity([profile_collect[ind] for ind in select_ind], mode="sum")
            TotInt_Collect.append(tot)
            # Draw profiles of FP distribution.
            #--> Calculate Polarity Score and cell length.
            #--> Process the profiles.
            EditProfiles = flip_profile(profile_collect, 1/3)
            EditProfiles = Batch_LenNorm(EditProfiles, option='points')
            avg_arr, std_arr, EditProfiles = profile_statistics(EditProfiles, Norm=False)
            #--> Single cell traces (marked one line).
            norm = max([max(arr) for arr in EditProfiles])
            self.PltFunc.ProfilePlot(EditProfiles, subfolder[0], mean_line=avg_arr, norm=norm)
            #--> Average of profiles with STD.
            self.PltFunc.ProfileStatPlot(avg_arr, std_arr, subfolder[0])

        # Draw statistics of total intensity.
        #--> Basic bar plot with mean and std.
        # Draw statistics of total intensity.
        #--> Basic bar plot with mean and std.
        max_bound = []
        for arr in TotInt_Collect:
            arr = arr[2]
            Q1 = np.percentile(arr, 75)
            Q3 = np.percentile(arr, 25)
            max_bound.append(Q1+1.5*(Q1-Q3))
        max_bound = max(max_bound)
        
        if 'tot-bar' in plots:
            self.PltFunc.BarPlot([v[0]/max_bound for v in TotInt_Collect],
                                 [v[1]/max_bound for v in TotInt_Collect],
                                 name_collect,
                                 'TotalIntensity',
                                 arrs=[v[2] for v in TotInt_Collect])
        if 'tot-box' in plots:
            #--> Box plot with quaters.
            self.PltFunc.BoxPlot([v[2]/max_bound for v in TotInt_Collect],
                                 name_collect,
                                 'TotalIntensity')
        if 'tot-violin' in plots:
            self.PltFunc.ViolinPlot([v[2]/max_bound for v in TotInt_Collect],
                                 name_collect,
                                 'TotalIntensity')
        
        df = pd.DataFrame(data=[v[2]/max_bound for v in TotInt_Collect], index=name_collect)
        df.to_csv('TotalIntensity.csv')
        print('The analysis is finishing.')


    def TwoChannelPipline(self, plots, twoselect=['', ''], FPextension=[-1, -1]):
        data = self.Data
        ChNum = self.ChNum
        # The essential containers.
        name_collect = [] # For storing the name of experiments.
        item = list(data.items()) # Create an iterable object from the input dictionary.
        TotInt_Collect = [] # For storing processed data related to total intensity.
        PS_collect = [] # For collecting data processed by Polarity Score Method.
        MC_collect = [] # For collecting data processed by Mass Center Method.
        corr_collect = [] # Confirm the colocalization by pearson correlation.
        arr1collect = []
        arr2collect = []
        LR_collect = []
        LR_max_collect = []
        Pearson_collect = []
        extension_saver = []
        for GroupNum, subfolder in enumerate(item): # Start the FOR LOOP by dictionary's KEYs.
            unipolar_index = [] # For collecting the cell ID whose PopZ is unipolar.
            CH1_index = [] # For collecting the location of PopZ in the cell.
            CH1_profiles = [] # For collecting the unipolar PopZ profiles.
            profile_collect = [] # For collecting the 2nd channel profiles with unipolar PopZ.
            PR_collect = [] # For collecting the index of peak of every single cells.
            for i in range(len(subfolder[1])-2):
                # Collect the PopZ or Channel1 raw profile.
                profile = subfolder[1]['{}'.format(i)]['S0']['profile']
                # Classify the spatial pattern of cells.
                if GroupNum == FPextension[1]:
                    if i in extension_saver:
                        polarity = ('unipolar', '')
                    else:
                        polarity = ('', '')
                else:
                    polarity = PolarClassifier(profile, 1/3)
                # Further analysis if unipolar
                if polarity[0] == 'unipolar' or subfolder[1]['PopZ']==0:
                    if GroupNum == FPextension[0]:
                        extension_saver.append(i)
                    # Record the peak index.
                    PR_collect.append(polarity[1])
                    # Record the unipolar cell index.
                    unipolar_index.append(i)
                    # Get and save the channel1 profiles.
                    CH1_profiles.append(subfolder[1]['{}'.format(i)]['S0']['profile'])
                    # Once there are more than 1 channel that are waiting for analysis...
                    if ChNum > 1:
                        try:
                            # Save the profile from getting channel2 data.
                            profile_collect.append(subfolder[1]['{}'.format(i)]['S1']['profile'])
                        except:
                            # When the .mat file contains some cells with ONLY ONE channal,
                            # the conditional command will skip the cell and print an error report.
                            print('[Data Error] The channel number does not match the setting (Set=2, D<2).',
                          subfolder[0], subfolder[1]['{}'.format(i)].keys())
                    # Find the biggest pole position in the cells.
                    if subfolder[1]['PopZ']==0:
                        pos, = np.random.choice([0,1], 1)
                        CH1_index.append(pos)
                    else:
                        if np.mean(profile[:int(len(profile)/2)])>np.mean(profile[int(len(profile)/2):]):
                            CH1_index.append(0) # 1: PopZ pole at left
                        else:
                            CH1_index.append(1) # 0: PopZ pole at right
            if len(CH1_index)>180:
                sampnum = np.random.choice(np.arange(130, 150))
                select_ind = np.random.choice(np.arange(len(CH1_index)), sampnum, replace=False)
            else:
                select_ind = np.arange(len(CH1_index))
            select_ind = select_ind.astype(int)
            print('The group', subfolder[0], 'has n=', len(select_ind), 'after unipolar PopZ selection.')
            with open('SampleNumber.txt', 'a+') as f:
                f.write('{0} n={1} \n'.format(subfolder[0], len(select_ind)))
            # Collect the name of the plasmid that transformed the cell.
            name_collect.append(subfolder[0])
            # Compute the statistics of total FP intensity.
            tot = TotalIntensity([profile_collect[ind] for ind in select_ind], mode="sum")
            TotInt_Collect.append(tot)
            # This part start to calculate the quantitative properties of the polar signal.
            #--> Compute the Polarity Score Method
            PS = PolarityScore([profile_collect[ind] for ind in select_ind], PopZ=CH1_index)
            PS_collect.append((PS[0], PS[1], PS[2], PS[3])) # 0: length, 1: single values, 2: mean, 3: std
            #--> Compute the mass center score of all profiles.
            MC = MassCenterScore([profile_collect[ind] for ind in select_ind], PopZ=CH1_index)
            MC_collect.append((MC[0], MC[1], MC[2], MC[3])) # 0: length, 1: single array, 2: mean, 3: std
            #--> Compute the difference between the signal in both poles.
            LR_collect.append(LR_diff([profile_collect[ind] for ind in select_ind], PopZ=CH1_index, ratio=1/3))
            # Pearson correlation coefficient.
            pc = Pearson_Correlation([profile_collect[ind] for ind in select_ind],
                                    [CH1_profiles[ind] for ind in select_ind])
            Pearson_collect.append(pc)
            # Process the profiles.
            #--> Process the PopZ profile.
            PProfiles = flip_profile(CH1_profiles, 1/3)
            PProfiles = Batch_LenNorm(PProfiles, option='points')
            #--> Process the 2nd channel profile.
            EditProfiles = flip_profile(profile_collect, 1/3, PopZ=CH1_index)
            EditProfiles = Batch_LenNorm(EditProfiles, option='points')
            EditProfiles_tmp = EditProfiles
            #--> Compute the statistic result of the profiles.
            _, _, PProfiles = profile_statistics(PProfiles, Norm=True, NormV=6.797564301714801)
            _, _, EditProfiles_tmp = profile_statistics(EditProfiles_tmp, Norm=False)
            
            if 'twoprofiles' in plots:
                # Save sample number.
                print('The group', subfolder[0], 'has n=', len(select_ind))
                with open('TwoProfiles_SampleNumber.txt', 'a+') as f:
                    f.write('{0} n={1} \n'.format(subfolder[0], len(select_ind)))
                # Save profiles for comparison.
                if subfolder[0]==twoselect[0]:
#                     print('save!!')
                    avg_arr, std_arr, EditProfiles = profile_statistics([
                        EditProfiles[ind] for ind in select_ind], Norm=False)
                    arr1collect = (avg_arr, std_arr)
                    
                if subfolder[0]==twoselect[1]:
#                     print('save!!')
                    avg_arr, std_arr, EditProfiles = profile_statistics([
                        EditProfiles[ind] for ind in select_ind], Norm=False)
                    arr2collect = (avg_arr, std_arr)
            else:
                avg_arr, std_arr, EditProfiles = profile_statistics(EditProfiles, Norm=False)
                    
            # Draw profiles of downstream FP distribution.
            if 'trace' in plots:
                #--> Single cell traces (marked one line).
                norm = max([max(arr) for arr in EditProfiles_tmp])
                self.PltFunc.ProfilePlot(EditProfiles_tmp, subfolder[0], mean_line=avg_arr, norm=norm)
            if 'statprofile' in plots:
                #--> Average of profiles with STD.
                self.PltFunc.ProfileStatPlot(avg_arr, std_arr, subfolder[0])
            
        if 'twoprofiles' in plots:
            max1 = np.max(arr1collect[0])
            max2 = np.max(arr2collect[0])
            maxmax = np.max([max1, max2])
            self.PltFunc.TwoProfilePlot(arr1collect, arr2collect, twoselect)
            df1 = pd.DataFrame(data = arr1collect[0]/maxmax, columns=['{}-mean'.format(twoselect[0])])
            df2 = pd.DataFrame(data = arr1collect[1]/maxmax, columns=['{}-std'.format(twoselect[0])])
            df3 = pd.DataFrame(data = arr2collect[0]/maxmax, columns=['{}-mean'.format(twoselect[1])])
            df4 = pd.DataFrame(data = arr2collect[1]/maxmax, columns=['{}-std'.format(twoselect[1])])
            df = pd.concat([df1, df2, df3, df4], axis=1, ignore_index=True)
            df.columns = ['{}-mean'.format(twoselect[0]),
                          '{}-std'.format(twoselect[0]),
                          '{}-mean'.format(twoselect[1]),
                          '{}-std'.format(twoselect[1])]
            df.to_csv('Linescan_comparison.csv')
        # Draw statistics of total intensity.
        #--> Basic bar plot with mean and std.
        if 'tot' in plots:
            max_bound = []
            for arr in TotInt_Collect:
                arr = arr[2]
                Q1 = np.percentile(arr, 75)
                Q3 = np.percentile(arr, 25)
                max_bound.append(Q1+1.5*(Q1-Q3))

            max_bound = max(max_bound)
            self.PltFunc.ViolinPlot([v[2]/max_bound for v in TotInt_Collect],
                                 name_collect,
                                 'TotalIntensity')
            df = pd.DataFrame(data=[v[2]/max_bound for v in TotInt_Collect], index=name_collect)
            df.to_csv('TotalIntensity.csv')

        # Draw polarity analysis.
        if 'ps' in plots:
            #--> Polarity Score with bar plot.
            self.PltFunc.ViolinPlot([v[1] for v in PS_collect],
                                 name_collect,
                                 'PolarityScore')
        if 'ms' in plots:
            #--> Mass Center Distribution with Quaters.
            self.PltFunc.HORBoxPlot([v[1] for v in MC_collect],
                                    name_collect)

            df = pd.DataFrame(data=[v[1] for v in MC_collect], index=name_collect)
            df.to_csv('MassCenter.csv')

        if 'lr' in plots:
            #--> LR ratio plot
            self.PltFunc.ViolinPlot([v[1] for v in LR_collect],
                                 name_collect,
                                'LRratio')

            df = pd.DataFrame(data=[v[1] for v in LR_collect], index=name_collect)
            df.to_csv('LRratio.csv')
        
        if 'pearson' in plots:
            # Draw bar plot with std to show the Pearson correlation coefficient.
            self.PltFunc.ViolinPlot([v for v in Pearson_collect],
                                             name_collect,
                                            'Pearson')
            df = pd.DataFrame(data=[v for v in Pearson_collect], index=name_collect)
            df.to_csv('Pearson.csv')
            
        print('The analysis is finishing.')




    def TimeLapsePipline(self):
        data = self.Data
        # The essential containers.
        name_collect = [] # For storing the name of experiments.
        item = list(data.items()) # Create an iterable object from the input dictionary.
        
        for subfolder in item: # Start the FOR LOOP by dictionary's KEYs.
            unipolar_index = [] # For collecting the cell ID whose PopZ is unipolar.
            PopZ_index = [] # For collecting the location of PopZ in the cell.
            PopZ_profiles = [] # For collecting the unipolar PopZ profiles.
            profile_arrs = [] # For collecting the 2nd channel profiles with unipolar PopZ.
            t_record = []
            kymo_arrs = []
            raw_save = []
            mother_tot = np.zeros((len(subfolder[1])-2, 11))
            daughter_tot = np.zeros((len(subfolder[1])-2, 11))
            for i in range(len(subfolder[1])-2):
                # Collect the PopZ raw profile.
                profile_collect = {'ch1':[[], []], 'ch2':[[], []]}
                kymo_collect = {'ch1':[], 'ch2':[]}
                dict_time = subfolder[1]['{}'.format(i)]
                ctrl = 1
                for t in range(len(dict_time)):
                    x = len(dict_time)-1
                    dict_lineage = dict_time['T{}'.format(t)]
                    for l in range(len(dict_lineage)):
                        if l==1 and ctrl==1:
                            t_record.append((t, len(dict_time)-1))
                            ctrl = 0
                            
                        if l==2:
                            if self.ChNum==1:
                                ch1_profile = np.concatenate((dict_lineage['L{}'.format(1)]['S0']['profile'], dict_lineage['L{}'.format(2)]['S0']['profile']))
                                ch2_profile = np.concatenate((dict_lineage['L{}'.format(1)]['S0']['profile'], dict_lineage['L{}'.format(2)]['S1']['profile']))
                            else:
                                ch1_profile = np.concatenate((dict_lineage['L{}'.format(1)]['S0']['profile'], dict_lineage['L{}'.format(2)]['S0']['profile']))
                                ch2_profile = np.concatenate((dict_lineage['L{}'.format(1)]['S1']['profile'], dict_lineage['L{}'.format(2)]['S1']['profile']))
                                
                            profile_collect['ch1'][1] = ch1_profile
                            profile_collect['ch2'][1] = ch2_profile

                        else:
                            if self.ChNum==1:
                                ch1_profile = dict_lineage['L{}'.format(l)]['S0']['profile'][0][0].flatten()
                                ch2_profile = dict_lineage['L{}'.format(l)]['S0']['profile'][0][0].flatten()
                            else:
                                ch1_profile = dict_lineage['L{}'.format(l)]['S0']['profile'][0][0].flatten()
                                ch2_profile = dict_lineage['L{}'.format(l)]['S1']['profile'][0][0].flatten()
                            profile_collect['ch1'][l] = ch1_profile
                            profile_collect['ch2'][l] = ch2_profile
                        
                    kymo_tmp1, kymo_tmp2, mother, daughter = kymo_preprocess(profile_collect['ch1'], profile_collect['ch2'])
                    daughter_tot[i, t] = daughter
                    mother_tot[i, t] = mother
                    kymo_collect['ch1'].append(kymo_tmp1)#/kymo_tmp1.max())
                    kymo_collect['ch2'].append(kymo_tmp2)#/kymo_tmp2.max())

                kymo_arrs.append(kymo_collect)
            mmax = mother_tot.max(axis=1)
            daughter_tot = np.array([d/mm for d, mm in zip(daughter_tot, mmax)])
            mother_tot = np.array([m/mm for m, mm in zip(mother_tot, mmax)])

            max_divide_t = np.max(np.array([t[0] for t in t_record]))
            max_final_t = np.max(np.array([t[1]-t[0] for t in t_record]))
            min_divide_t = np.min(np.array([t[0] for t in t_record]))
            min_final_t = np.min(np.array([t[1]-t[0] for t in t_record]))
            
            for t_points, kymodict in zip(t_record, kymo_arrs):
                aaa=kymodict['ch1']
                kymodict['ch1'] = [np.array([]) for _ in range(max_divide_t-t_points[0])]+kymodict['ch1']+[np.array([]) for _ in range(max_final_t-(t_points[1]-t_points[0]))]
                kymodict['ch2'] = [np.array([]) for _ in range(max_divide_t-t_points[0])]+kymodict['ch2']+[np.array([]) for _ in range(max_final_t-(t_points[1]-t_points[0]))]

            lapse_graph_ch1 = []
            lapse_graph_ch2 = []
            global_maxLen = []
            for t in range(len(kymo_arrs[0]['ch1'])):
                try:
                    [len(kymodict['ch1'][t]) for kymodict in kymo_arrs]
                except:
                    print(len(kymo_arrs[0]['ch1']), [len(kymodict['ch1']) for kymodict in kymo_arrs])
                max_len = np.mean(np.array([len(kymodict['ch1'][t]) for kymodict in kymo_arrs]))
                EdditkymoCH1 = []
                EdditkymoCH2 = []
                for kymodict in kymo_arrs:
                    if np.any(kymodict['ch1'][t]):
                        EdditkymoCH1.append(LenNorm(kymodict['ch1'][t], Len=int(2*max_len))[0])
                        EdditkymoCH2.append(LenNorm(kymodict['ch2'][t], Len=int(2*max_len))[0])
                lapse_graph_ch1.append(np.mean(EdditkymoCH1, axis=0))
                lapse_graph_ch2.append(np.mean(EdditkymoCH2, axis=0))
                global_maxLen.append(int(2*max_len))
            global_maxLen = np.max(np.array(global_maxLen))
            Final_ch1 = np.zeros((len(lapse_graph_ch1), global_maxLen))*0
            Final_ch2 = np.zeros((len(lapse_graph_ch1), global_maxLen))*0
            bound_x = np.array([np.zeros(len(lapse_graph_ch1)), np.zeros(len(lapse_graph_ch1))])
            bound_y = np.zeros(len(lapse_graph_ch1))
            for ind in range(len(lapse_graph_ch1)):
                boundary = int((len(Final_ch1[0, :])-len(lapse_graph_ch1[ind]))/2)
                Final_ch1[ind, boundary:boundary+len(lapse_graph_ch1[ind])] = lapse_graph_ch1[ind]
                Final_ch2[ind, boundary:boundary+len(lapse_graph_ch2[ind])] = lapse_graph_ch2[ind]
                bound_x[0][ind] = boundary
                bound_x[1][ind] = len(Final_ch1[0, :])-boundary-1
                bound_y[ind] = ind
            bound = (bound_x, bound_y)
            daughter_tot[:, 0:5] = mother_tot[:, 0:5]
            
            Final_ch2 = Final_ch2[:10, :]
            print('The maximum value is', np.max(Final_ch2))
            Final_ch2 = Final_ch2/np.max(Final_ch2)
            self.PltFunc.divergent(mother_tot, daughter_tot)
            self.PltFunc.kymograph(Final_ch1, max_divide_t, bound)
            df = pd.DataFrame(Final_ch2)
            df.to_csv('kymograph_raw.csv')
            self.PltFunc.kymograph(Final_ch2, max_divide_t, bound)
            
            ind1 = 0
            ind2 = 0
            data = np.zeros((daughter_tot.shape[0]*2, daughter_tot.shape[1]))
            for i in range(daughter_tot.shape[0]*2):
                if i%2==0:
                    data[i, :] = daughter_tot[ind1, :]
                    ind1+=1
                else:
                    data[i, :] = mother_tot[ind2, :]
                    ind2+=1
            
            ID = np.arange(daughter_tot.shape[0])
            sub_col = ['daughter', 'mother']
            iterables = [ID, sub_col]
            arrays = pd.MultiIndex.from_product(iterables, names=['ID', 'Cell type'])
            df = pd.DataFrame(data, index=arrays)
            df = df.T
            df.to_csv('timelapse_raw.csv')
            print('The analysis is finishing.')