import pandas as pd
import os


homedir='/Users/owner/Functional_Connectivity'
coordfile='Template.nii'
radius=25

df=pd.read_csv(homedir+'/'+coordfile) #read the csv and save it as a pandas dataframe

#loop through each row of the dataframe, create variables for x,y,z @ each row, run 3dcalc inserting these variable
for row in df.iterrows():
    x_coord=row[1][0]*-1
    y_coord=row[1][1]*-1
    z_coord=row[1][2]
    outname=row[1][3]
    
    command="""3dcalc -a 1mm_ch2better-2.nii -expr 'step(%s-(x-%s)*(x-%s)-(y-%s)*(y-%s)-(z-%s)*(z-%s))' -prefix %s""" %(radius, x_coord, x_coord, y_coord, y_coord, z_coord, z_coord, outname+'.nii')
    os.system(command)


#thats it so far