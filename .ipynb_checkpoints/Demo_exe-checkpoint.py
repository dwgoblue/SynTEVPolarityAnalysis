import time
from Piplines import*


def Main():
    # Apply absolute path if needed.
    path = r'./Demo/[Demo]DataInput_folder.txt'

    # Change the first input from _ to "" in older version.
    Starter = PiplinesObj("", TimeLapse=False,
                          columns=[(0,1)],
                          mode='txt',
                          txt_path=path)

    # Analysis of LR ratio and two-profiles comparison.
    Starter.TwoChannelPipline(plots=['lr', 'twoprofiles'],
                              twoselect=['NSSCP-32', 'NSSC-32'])
    
    time.sleep(5)
    
if __name__ == '__main__':
    Main()
    
    
    