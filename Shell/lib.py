import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
# from scipy.optimize import curve_fit
from pathlib import Path
from math import sqrt,sinh,asinh,exp,acos,pi,asin,tanh

plt.style.use('classic')

mypara = {'axes.labelsize': 30,
            'axes.titlesize': 15,
            'axes.facecolor': 'white',
            # 'font.size': 14,
            'legend.fontsize': 15,
            'xtick.labelsize': 29,
            'ytick.labelsize': 29,
            'font.family': 'times new roman',
            'figure.figsize': (8, 6)}

plt.rcParams.update(mypara)
fig, ax = plt.subplots()   # maybe need to be modified
def fig_decor(title,xlab,ylab):
    plt.xlabel(xlab)
    plt.ylabel(ylab)
    plt.title(title)
    plt.legend(loc='upper right', frameon=False)
    plt.minorticks_on()
    # plt.yscale('linear')
    # plt.ticklabel_format(useOffset=False, style='plain')
    plt.grid(False)
    plt.tight_layout()

colors = cm.tab20(np.linspace(0, 1, 20))
def raw_data(data, yerr=None, title='', xlabel='', ylabel='', labels=None,filep='comp.png'):
    for i, (x, y) in enumerate(data):
        label = labels[i] if labels is not None else f'Line {i}'
        color = colors[i % 20]  
        if yerr is not None:
            ax.errorbar(x, y, yerr=yerr[i], label=label, fmt='o-', color=color)
        else:
            ax.plot(x, y, label=label)  # here need to be modified
    fig_decor(title,xlabel,ylabel)
    plt.savefig(filep, format='png', dpi=900)
    # plt.show()

def data_process(index,index2,filen='output.dat',hqg='cl_temp'):
    # bootstarp function
    nm=Path(filen).stem
    fileparent=Path('processed_data/'+hqg)
    if not fileparent.exists():
        fileparent.mkdir()
    file=fileparent/(nm+str(index)+'.txt')
    if file.exists():
        print('file exists')
        return file
    filepath = sorted(Path('.').rglob(filen))
    data0list=[]
    for i in filepath:
        print(i)
        data=np.loadtxt(i,usecols=(0,index2,index)) #revise here
        data0samples=[]
        bins=len(data)
        data0samples=data 
        # if filen=='out.dat' and index==6:data0samples[:,2]=(2-data[:,2])
        # if filen=='out.dat' and index==10:data0samples[:,2]=2-data[:,2]
        # if filen=='out.dat' and index==12:data0samples[:,2]=(5-3*data[:,2])/2
        if data.ndim==2:
            data0list.append([data[0][0],data[0][1],data0samples.mean(0)[2],sqrt(data0samples.var(0)[2]/(bins-1))])
        else:
            data0list.append([data[0],data[1],data[2],0])
    np.savetxt(file,sorted(data0list,key=lambda x:x[0]),fmt='%d'+'   %.5f'+2*'   %.12f')
    return file

def plot_data(file_name, file_saved, ele_values, istemp_index):
    file_path = data_process(istemp_index[0],istemp_index[1], filen=file_name, hqg=file_saved)
    data_array = np.loadtxt(file_path)
    
    ele_list = [data_array[np.where(data_array == ele)[0]] for ele in ele_values]
    data, yerr, labels = [], [], []
    
    for i, ele in enumerate(ele_list):
        data.append((ele[:, istemp_index[2]], ele[:, 2]))
        yerr.append(ele[:, 3])
        labels.append(f'$N$={ele_values[i]}')
    
    raw_data(data, yerr=yerr, title='', xlabel='$?$', ylabel=r'$?$', labels=labels, filep=file_path.with_suffix('.png'))

from PIL import Image

def convert_png_to_eps(png_file, eps_file):
    with Image.open(png_file) as im:
        im.save(eps_file, 'eps')

def rebin(re,filen):
    filepath = sorted(Path('.').rglob(filen))
    for i in filepath:
        data=np.loadtxt(i)
        if len(data)>re:
            ans=[data[i*re:(i+1)*re,:].mean(0) for i in range(len(data)//re)]
            np.savetxt(i,ans,fmt='%d'+13*'   %.12f',delimiter='   ')
        else:
            print(i)

def binder(pos,n,tc,v):
    dt = np.loadtxt('processed_data/'+pos+'/out4.txt')
    dt[:,1] = (dt[:,1] - tc) * dt[:,0]**(1/v)
    for i in n:
        ele = dt[np.where(dt[:,0]==i)[0]]
        plt.errorbar(ele[:,1],ele[:,2],ele[:,3],fmt='o',linestyle='-',label='N={}'.format(i),markersize=3)
    plt.legend(loc ='best', frameon=False)
    plt.show()

def m_2(pos,n,tc,v,em):
    dt = np.loadtxt('processed_data/'+pos+'/out6.txt')
    dt[:,2] = dt[:,2] * dt[:,0]**(em) #- 4.3*dt[:,0]**(-0.91/2) #+ 2.9*(dt[:,1] - 0.947)*dt[:,0]**(-0.05)
    dt[:,1] = (dt[:,1] - tc) * dt[:,0]**(1/v)
    for i in n:
        ele = dt[np.where(dt[:,0]==i)[0]]
        plt.plot(ele[:,1],ele[:,2],label='N={}'.format(i),marker='o',markersize=3)
    # plt.legend()
    plt.show()
