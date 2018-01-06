#!/usr/bin/python

import sys
import os
import string
import csv
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import ipc
from bioservices import UniProt

def IndexMaker(x,y=0,way='l'):
    index=[]
    if way == 'l':
        for l in list(string.ascii_uppercase[:x]):
            index.append(l)
        return index
    if way =='ln':
        for l in list(string.ascii_uppercase[:x]):
            for n in range(1,y+1):
                index.append(l+str(n))
        return index
    else:
        print("Error in function IndexMaker!! please define if you want numbers (IndexMaker(x,'l')) or letternumber combination (IndexMaker(x,'ln')) ")
        sys.exit(1)

def CSV_cutter(csvfile,cut_folder,sig_prot):
    filename=csvfile.rsplit('/',1)[1]
    with open(csvfile,'r') as csvin:
        csvreader= csv.reader(csvin,delimiter=',')
        Writing=False
        first_iteration=True
        #Following loop goes though csv file
        #start writing when array parts ir reached
        #Stop when all the useful stuff is read
        for line in csvreader:
            if not line: #Skips empty lines
                continue
            if line[0:1]==['prot_hit_num']: #On-switch
                Writing=True
                width=len(line)+2
            if line[0:1]==['Queries']: #Off-switch
                break
            if Writing == True:
                if first_iteration==True: #This writes heades in so it fits with the emPAI row
                    cont=line+['emPAI','emPAI_value']
                    columns=cont
                    df=pd.DataFrame(columns=columns)
                    first_iteration=False
                elif len(line) < width:
                    cont=line+['NaN','NaN']
                    df=df.append(pd.DataFrame(cont,index=columns).T,ignore_index=True)
                elif len(line) >=width:
                    df=df.append(pd.DataFrame(line,index=columns).T,ignore_index=True)
            else:
                continue
    df_colcut=df.loc[:,['prot_acc','prot_desc','prot_score','prot_mass','prot_sequences_sig','emPAI_value']]
    #df_colcut['prot_acc']=df_colcut['prot_acc'].map(lambda x: x.lstrip('1::'))
    #df_colcut['prot_acc']=df_colcut['prot_acc'].map(lambda x: x.lstrip('2::')) 
    #df_colcut['prot_acc']=df_colcut['prot_acc'].map(lambda x: x.lstrip('3::'))
    df_colcut['emPAI_value']=df_colcut['emPAI_value'].apply(pd.to_numeric, errors='coerce')
    df_colcut['prot_sequences_sig']=df_colcut['prot_sequences_sig'].apply(pd.to_numeric, errors='coerce')
    df_output=df_colcut.loc[(df_colcut['emPAI_value'] > 0) & (df_colcut['prot_sequences_sig'] > sig_prot)]
    df_output.to_csv(cut_folder+'/cut_%s'% filename,index=False)
    return df_output


def Categorizer(dataframe,x,y):
    df=dataframe
    lowerlimit=int(x*0.8)
    DFcoord=df[df.columns[df.columns.get_loc("A1"):]]
    TotalPixl=x*y
    upperlimit=int(round(TotalPixl/2))
    Summer=TotalPixl-DFcoord[DFcoord > 0].isnull().sum(axis=1)
    df['num_pixels']=Summer
    df.loc[df['num_pixels'] < lowerlimit, 'Category']="found"
    df.loc[df['num_pixels'] > upperlimit, 'Category']="overhalf"
    df.loc[(df['num_pixels'] >= lowerlimit)&(df['num_pixels']<=upperlimit), 'Category']="useful"
    return df


def RowNormalizer(dataframe):
    df=dataframe
    DFcoord=df[df.columns[df.columns.get_loc("A1"):]]
    bf=df[df.columns[:df.columns.get_loc("A1")]]
    normDFcoord=DFcoord.div(DFcoord.max(axis=1), axis=0)
    df=pd.concat([bf,normDFcoord], axis=1)
    return df



def CoordinateExtractor(dataframe,x,y):
    df=dataframe
    letterindex=IndexMaker(x,'l')
    numbercols=np.arange(1,y+1)
    df=df.astype('object')
    for index,row in df.iterrows():
        rsdf=pd.DataFrame(row.loc['A1':].values.reshape((x,y)),
                          index=letterindex, columns=numbercols, dtype='float64').T
        rsdf_V_flip=np.flipud(rsdf.values)
        nps=np.argwhere(rsdf_V_flip >0)
        coord_int=[]
        for i in nps:
            lst=[i[1]*(8/x)+((8/x)/2),i[0]]
            lst.append(rsdf_V_flip[i[0]][i[1]])
            coord_int.append(lst)
        df.loc[index,"Coordinates_emPAI"]=coord_int
    return df

def Smear_check(dataframe,y):
    df=dataframe
    df_t=df[df.columns[df.columns.get_loc("A1"):df.columns.get_loc("B1")]]
    null_count=df_t[df_t > 0].isnull().sum(axis=1)
    df["Smear"]=null_count < (y/2)
    return df

def Smear_check_th(dataframe,y):
    df=dataframe
    df_t=df["Threshold"]
    df_s=df[df.columns[df.columns.get_loc("A1"):df.columns.get_loc("B1")]]
    for index,row in df_s.iterrows():
        null_count=row.where(row > df_t[index]).isnull().sum()
        df.loc[index,"Smear"]=null_count < (y/2)
    return df

def pI_calc(dataframe):
    df=dataframe
    u = UniProt()
    for index,row in df.iterrows():
        seqce=u.search(df.loc[index,"prot_acc"],frmt="tab",columns="sequence").split('\n')
        p_i=ipc.predict_isoelectric_point(seqce[1])
        df.loc[index,"pI"]=p_i

    return df


def CSVconverter(input_folder,cut_folder,dbfile,x,y,taxonomy='',sig_seq=1):

    if not os.path.isdir(cut_folder):
        os.makedirs(cut_folder)

    real_list=[]

    comp_list=IndexMaker(x,y,'ln')
    prot_cols=['prot_acc','prot_desc','prot_mass','pI','[Urea]_50','Delta_G','Category',
               'num_pixels','Smear','Coordinates_emPAI','Threshold']
    columnsdb=prot_cols+comp_list

    df_db=pd.DataFrame(columns=columnsdb)
    df_db["Threshold"]= 0.0

    for entry in os.scandir(os.getcwd()+'/'+input_folder):
        if entry.name.startswith('.'):
            continue
        filename=entry.name
        
        filename=filename[:-4]
        real_list.append(filename)
        df_csv=CSV_cutter(entry.path,cut_folder,sig_seq)
        df_csv=df_csv[['prot_acc', 'prot_desc', 'prot_mass','emPAI_value']]
        df_csv.rename(columns = {'emPAI_value':filename}, inplace=True)
        df_db=pd.concat([df_db,df_csv],ignore_index=True).reset_index(drop=True)

    df_db=df_db[columnsdb].fillna(0)
    df=df_db.groupby(prot_cols).sum()
    df=df.reset_index()
    df['prot_desc'].replace(regex=True,to_replace='/',value='-',inplace=True)
    if taxonomy !='':
        taxstring="OS="+ taxonomy
        df=df[df["prot_desc"].str.contains(taxstring)]
        df.reset_index(inplace=True)
    df=Categorizer(df,x,y)
    df=Smear_check(df,y)
    df=CoordinateExtractor(df,x,y)
    df=pI_calc(df)
    df.to_csv(dbfile,index=False)


    differ_list=list(set(comp_list)-set(real_list))
    if differ_list:
        readfilename=os.getcwd()+'/'+cut_folder+"/cut_"+real_list[0]+".csv"
        df_miss=pd.read_csv(readfilename)
        df_empty=pd.DataFrame(columns=df.columns)
        for coord in differ_list:
            filename=coord+".csv"
            df_empty.to_csv(cut_folder+'/cut_%s'% filename,index=False)
            
         
    df_plot= RowNormalizer(df)
    DFcoord=df[df_plot.columns[df_plot.columns.get_loc("A1"):]]

    count_inverse_sum=DFcoord[DFcoord > 0].isnull().sum(axis=0)
    count_sum= DFcoord.index.size - count_inverse_sum

    count_rsdata=count_sum.values.reshape((x,y)).T

    letterindex=IndexMaker(x,'l')                                                              
    numbercols=np.arange(1,y+1)                                                                
    
    sns.set(style="white")
    sns.set(font_scale=1.3)
    fig, ax1 = plt.subplots(figsize=(12,9))


    df_count=pd.DataFrame(count_rsdata,columns=letterindex,index=numbercols)
 
    sns_plot=sns.heatmap(df_count,cmap="Blues",annot=True,cbar=True, fmt='d', ax=ax1)

    sns_plot.xaxis.tick_top()

    ax1.set_title("Protein hits. pr pixel",size=50,y=1.1)
    plt.tight_layout()
    plt.savefig("%s.png"%dbfile.split('_db')[0],dpi=500)
    plt.close(fig)        

    return df
            
def TrypsinCSVconverter(input_folder,cut_folder,dbfile,x,y,sig_seq=0):

    if not os.path.isdir(cut_folder):
        os.makedirs(cut_folder)

    comp_list=IndexMaker(x,y,'ln')
    letterindex=IndexMaker(x,'l')   
    prot_cols=['prot_acc','prot_desc','prot_score']
    columnsdb=prot_cols+comp_list

    df_db=pd.DataFrame(columns=columnsdb)
    

    names=[]
    
    for entry in os.scandir(os.getcwd()+'/'+input_folder):
        if entry.name.startswith('.'):
            continue
        filename=entry.name
        filename=filename[:-4]
        df_csv=CSV_cutter(entry.path,cut_folder,sig_seq)
        df_csv=df_csv[df_csv['prot_desc'].str.contains('TRYP_PIG',case=False)]
        if df_csv.empty:
            df_csv=pd.DataFrame([['P00761','TRYP_PIG Trypsin OS=Sus scrofa PE=1 SV=1',0,0]],
                             columns=['prot_acc', 'prot_desc','prot_score','emPAI_value'])
        else:
            df_csv=df_csv[['prot_acc', 'prot_desc','prot_score','emPAI_value']]
        df_csv.rename(columns = {'emPAI_value':filename}, inplace=True)
        df_db=pd.concat([df_db,df_csv],ignore_index=True).reset_index(drop=True)
        names.append(filename)
    
    
    df_db=df_db[columnsdb].fillna(0)
    df_db=df_db.reindex_axis(columnsdb,axis=1) 
    df_db.drop(['prot_score'],axis=1, inplace=True)
    df=df_db.groupby(['prot_acc','prot_desc']).sum()
    df=df.reset_index()

    
    numbercols=np.arange(1,y+1)                                                                
    sns.set(style="white") 
    sns.set(font_scale=1.3)

    Tryp_coord=df[df.columns[df.columns.get_loc("A1"):]]
    tryp_rsdata=Tryp_coord.values.reshape((x,y)).T

    df_tryp_rs=pd.DataFrame(tryp_rsdata,columns=letterindex,index=numbercols)
    
    fig, ax = plt.subplots(figsize=(12,9))
    sns_plot=sns.heatmap(df_tryp_rs,cmap="Blues",annot=True,cbar=True,fmt='.3g')
    sns_plot.xaxis.tick_top()
    
    ax.set_title("Trypsin emPAI pr. pixel",size=50, y=1.1)
    plt.tight_layout()
    plt.savefig("%s_emPAI.png"%dbfile.split('_db')[0],dpi=500)
    plt.close()
    
    
    return df