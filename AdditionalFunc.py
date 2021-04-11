import numpy as np
# Geometric methods
def centroid(vertexes):
    
    x = np.array([coord[0] for coord in vertexes])
    y = np.array([coord[1] for coord in vertexes])
    cx, cy, A = 0, 0, 0
    for i in range(len(x)-1):
        cx += (x[i] + x[i+1])*(x[i]*y[i+1]-x[i+1]*y[i])
        cy += (y[i] + y[i+1])*(x[i]*y[i+1]-x[i+1]*y[i])
        A += (x[i]*y[i+1]-x[i+1]*y[i])
    
    A = 0.5*A
    cx = cx/(6*A)
    cy = cy/(6*A)
    
    return np.array([cx, cy])

def linear_distance(a, b):
    return np.sqrt((a[0]-b[0])**2+(a[1]-b[1])**2)


def MassCenter1D(arr):
    mass = 0
    for i, v in enumerate(arr):
        mass += v*i
    if np.sum(arr)==0:
        return 0
    else:
        return mass/np.sum(arr)


def PeakIsolation(arr, max_ind, buffer):
    # Scan toward left.
    progress = True
    ind = 0
    tmp = arr[max_ind]
    left_mini = 0
    while progress==True:
        ind+=1
        
        if int(max_ind-ind)==-1:
            left_mini = int(max_ind-ind+1)
            #print('minimum at boundary.')
            left_len = ind-1
            progress = False
        
        elif arr[int(max_ind-ind)]>arr.mean() or arr[int(max_ind-ind)]<tmp:
            tmp = arr[int(max_ind-ind)]

        else:
            try:
                B = arr[int(max_ind-ind-np.arange(buffer))]
            except:
                B = arr[:int(max_ind-ind)]
            if np.all(B>=tmp):
                left_mini = int(max_ind-ind+1)
                left_len = ind-1
                progress = False

    # Scan toward right.
    progress = True
    ind = 0
    tmp = arr[max_ind]
    right_mini = 0
    while progress==True:
        ind+=1
        if int(max_ind+ind)==len(arr)-1:
            #print('right minimum at boundary!!!')
            right_mini = int(max_ind+left_len)
            progress = False
        
        elif arr[int(max_ind+ind)]<tmp or arr[int(max_ind+ind)]>arr.mean():
            tmp = arr[int(max_ind+ind)]

        else:
            try:
                B = arr[int(max_ind+ind+np.arange(buffer))]
                #print(B)
            except:
                B = arr[int(max_ind+ind):]
            if np.all(B>=tmp):
                if ind-1 < left_len:
                    right_mini = int(max_ind+left_len)
                else:
                    right_mini = int(max_ind+ind-1)
                progress = False

    return (left_mini, right_mini)       
    #return (int(max_ind-5), int(max_ind+5))
    
    
    
def confidence_ellipse(x, y, ax, n_std=3.0, facecolor='none', **kwargs):
    """
    Create a plot of the covariance confidence ellipse of `x` and `y`

    Parameters
    ----------
    x, y : array_like, shape (n, )
        Input data.

    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.

    n_std : float
        The number of standard deviations to determine the ellipse's radiuses.

    Returns
    -------
    matplotlib.patches.Ellipse

    Other parameters
    ----------------
    kwargs : `~matplotlib.patches.Patch` properties
    """
    if len(x) != len(y):
        raise ValueError("x and y must be the same size")

    cov = np.cov(x, y)
    pearson = cov[0, 1]/np.sqrt(cov[0, 0] * cov[1, 1])
    # Using a special case to obtain the eigenvalues of this
    # two-dimensionl dataset.
    ell_radius_x = np.sqrt(1 + pearson)
    ell_radius_y = np.sqrt(1 - pearson)
    ellipse = Ellipse((0, 0),
        width=ell_radius_x * 2,
        height=ell_radius_y * 2,
        facecolor=facecolor,
        **kwargs)

    # Calculating the stdandard deviation of x from
    # the squareroot of the variance and multiplying
    # with the given number of standard deviations.
    scale_x = np.sqrt(cov[0, 0]) * n_std
    mean_x = np.mean(x)

    # calculating the stdandard deviation of y ...
    scale_y = np.sqrt(cov[1, 1]) * n_std
    mean_y = np.mean(y)

    transf = transforms.Affine2D() \
        .rotate_deg(45) \
        .scale(scale_x, scale_y) \
        .translate(mean_x, mean_y)

    ellipse.set_transform(transf + ax.transData)
    return ax.add_patch(ellipse)
