
# coding: utf-8

# In[ ]:


#cd ~/Documents/molecules/p8git/results
get_ipython().magic(u'matplotlib inline')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
from scipy import interpolate
import IPython
import matplotlib.cm as cm

#!ls
def solid_read(file):
    return pd.read_table(file,delimiter=' ',index_col=0,names=['x','y','z','sxb','syb','szb','sxk','syk','szk','temp']).dropna(axis=1).reset_index()

def liquid_read(file):
    return pd.read_table(file,delimiter=' ',index_col=0,names=['id', 'x', 'y', 'z', 'Sx', 'Sy', 'Sz', 'Sxb', 'Syb', 'Szb', 'Sxk', 'Syk', 'Szk', 'Sxy','Syz','Sxz','Sxyb', 'Syzb','Sxzb', 'Sxyk','Syzk', 'Sxzk' ]).dropna(axis=1).reset_index()
#Read files -- reference and composite 

# K and epsilon values to be read
K_values=np.array(['5','10','107'])
e_values=np.array(['05','14','22'])
refs=np.array([]);
systems=np.array([]);
for k in K_values:
    str1="results/ref/stress.{}.out"
    f=str1.format(k)
    var1="ref_{}"
    var1=var1.format(k);
    refs=np.append(refs,var1)
    exec(var1 + "=solid_read(f)")
    for e in e_values:
        str2="results/composite/stress.{}.{}.out"
        ff=str2.format(k,e)
        var2="lq_{}_{}"
        var2=var2.format(k,e)
        exec(var2 + "=liquid_read(ff)")
        systems=np.append(systems,var2)


# In[ ]:


## Stick to scatter plots for now 
## Plot layers, bond/pair/ke, anisotropy --- at one value of K and epsilon first
## Find ways to quantify the decay function 

#Trail Plot
## Define contour plot custom made
## variables are input (lq or ref), quantity(sxb,syb), which layers to plotted, frames

def mycont(data,quant='sxb',layers=np.array([10,1]),nf=1):
    yr=layers.size*2
    #print layers.size
    fig,axs=plt.subplots(layers.size,1,figsize=(28,yr))
    j=0;
    for i in layers:
        l=np.arange((i-1)*3600,i*3600);
        
        f=axs[j].scatter(data.loc[l]['x'],data.loc[l]['y'],c=data.loc[l][quant],s=64,marker='h',cmap="jet",alpha=0.5)
        mean_z=data.loc[l]['z'].mean()
        axs[j].set_title(str(mean_z))
        axs[j].set_xlim([0,375])
        axs[j].set_ylim([0,60])
        j=j+1;
    cbar = fig.colorbar(f, ax=axs.ravel().tolist(), shrink=0.95)
#mycont(ref_10,quant='szb')  


# In[ ]:


#based on frames
def mycont_frame(data,quant='sxb',nl=1,nf=1,yr=6):
    #print layers.size
    fig,axs=plt.subplots(nf,1,figsize=(28,yr), squeeze=False)
    axs = axs.flatten()
    #print axs
    j=0;
    key=np.arange(nf);
    key=np.repeat(key,36000);
    trj=data.groupby(key);
    
    for i in np.arange(nf):
        
        #l=np.arange((nl-1)*3600,nl*3600);
        
        data_temp=data.loc[trj.groups[i]];
        keys=np.arange(10);
        keys=np.repeat(keys,3600);
        data_layers=data_temp.groupby(keys);
        
        data_temp=data_temp.loc[data_layers.groups[nl]];
        
        #print data_temp
        f=axs[j].scatter(data_temp['x'],data_temp['y'],c=data_temp[quant],marker='h',cmap='jet',s=64)
        mean_z=data_temp['z'].mean()
        axs[j].set_title(str(mean_z))
        axs[j].set_xlim([0,375])
        axs[j].set_ylim([0,60])
        j=j+1;
    cbar = fig.colorbar(f, ax=axs.ravel().tolist(), shrink=0.95)
    #plt.colorbar(f)
#lq_5_14.loc[l1]['sxb']## Read lammpstrj file 
#%automagic
#trj=pd.read_table('stress_conv/stress_new_all.out',delimiter=' ')
def solid_lmpread(file,nf=4,nsolid=36000):
    skiprows=np.array([],dtype=int);
    for i in range(nf):
        j=i*(nsolid+9);
        skiprows=np.append(skiprows,np.arange(j,j+9))
        #print skiprows
    names=np.array(['id', 'x', 'y', 'z','Sxb', 'Syb', 'Szb', 'Sxk', 'Syk', 'Szk']);
    return pd.read_table(file,delimiter=' ',skiprows=skiprows,names=names,index_col=False).dropna(axis=1)
    #key=np.arange(nf)
    #key=np.repeat(key,nsolid)
    #frames=trj.groupby(key)
#trj.loc[frames.groups[0]]['sxb']

def liquid_lmpread(file,nf=10,nsolid=36000):
    skiprows=np.array([],dtype=int);
    for i in range(nf):
        j=i*(nsolid+9);
        skiprows=np.append(skiprows,np.arange(j,j+9))
        #print skiprows
    names=np.array(['id', 'x', 'y', 'z', 'Sx', 'Sy', 'Sz', 'Sxb', 'Syb', 'Szb', 'Sxk', 'Syk', 'Szk', 'Sxy','Syz','Sxz','Sxyb', 'Syzb','Sxzb', 'Sxyk','Syzk', 'Sxzk']);
    return pd.read_table(file,delimiter=' ',skiprows=skiprows,names=names,index_col=False).dropna(axis=1)
    #key=np.arange(nf)
    #key=np.repeat(key,nsolid)
    #frames=trj.groupby(key)
#trj.loc[frames.groups[0]]['sxb']


# In[ ]:


## Frame averages

## nvt liquid
def frame_avg(data,nf,natom=36000):
    frames=np.arange(nf);
    keys=np.repeat(frames,natom)

    grouped=data.groupby(keys)
    dfs={}
    for i in range(nf): 
        dfs[i]=data.loc[grouped.groups[i]].reset_index(drop=True)

    panel=pd.Panel(dfs)
    return panel.mean(axis=0)

nvt_avg=frame_avg(nvt,10)
nvt_avg=frame_avg(npt,10)
nvt_solid_avg=frame_avg(solid_nvt,4)
npt_solid_avg=frame_avg(solid_npt,4)

