import os


#/data/soft/mricron/templates	 
atlas='/Users/rsblumenfeld/work/Atlases/Crossley/template.nii.gz' #change to reflect your path
#atlas='/soft/fsl/data/atlases/HarvardOxford/HarvardOxford-cort-maxprob-thr0-1mm.nii.gz'

#regions =[0:47]
labels = range(1,639)


for x,roi in enumerate(labels):
     out='C'+roi+'.nii.gz'	 	 
     afni_expression="3dcalc -a %s -expr 'equals(a,%s)' -prefix %s" %(atlas, int(roi), out)                                                
     #print afni_expression
     os.system(afni_expression)
