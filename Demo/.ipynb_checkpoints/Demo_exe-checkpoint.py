import time

def Main():
    from Piplines import*
    path = r'/home/r04b43015/SynPolarityAnalysis/Demo/[Demo]DataInput_folder.txt'
    Starter = PiplinesObj(_, TimeLapse=False,
                          columns=[(0,1)],
                          mode='txt',
                          txt_path=path)
    Starter.TwoChannelPipline(plots=['lr', 'twoprofiles'],
                              twoselect=['NSSCP-32', 'NSSC-32'])
    
    time.sleep(5)
    
if __name__ == '__main__':
    Main()
    
    
    