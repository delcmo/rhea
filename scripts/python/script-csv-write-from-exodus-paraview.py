#### import the simple module from the paraview
from paraview.simple import *
#### disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

import os, glob
# set path of the working directory
path_in='/Users/mxd/Output_files/rhea/mach-3-pure-abs-const-xs-ev' # path to be specified by users
os.chdir(path_in)
# set path of the target directory TODO
#path_out='' # path to be specified by users
print '#### Check current path: ####'
print os.path.realpath('.')
# set variables for loop
nel_start=200
nel_end=2100
increment=100

# loop over exodus files
print '#### Starting loop over exodus files ####'
for x in range(nel_start, nel_end+increment, increment):
  input_filename="mach-3-nb-cells-" + str(x) + "_out.e" # filenames should be changed to match user's needs
  output_filename="mach-3-nel-" + str(x) + "-points.csv" # filenames should be changed to match user's needs
  # create a new 'ExodusIIReader'
  reader = ExodusIIReader(FileName=input_filename)
  # get animation scene
  animationScene1 = GetAnimationScene()
  # update animation scene based on data timesteps
  animationScene1.UpdateAnimationUsingDataTimeSteps()
  # select all variables
  reader.SelectAllVariables()
  # go to last time step
  animationScene1.GoToLast()
  # save data
  SaveData(os.path.join(path_in, output_filename), proxy=reader, Precision=20)

print '#### Ending loop over exodus files ####'

# remove all files finishing with '1.csv', '2.csv', '3.csv' and '4.csv'
print '#### Removing all unecessary files (*1.csv, *2.csv, *3.csv and *4.csv) ####'
file1csv=glob.glob(os.path.join(path_in, '*1.csv'))
file2csv=glob.glob(os.path.join(path_in, '*2.csv'))
file3csv=glob.glob(os.path.join(path_in, '*3.csv'))
file4csv=glob.glob(os.path.join(path_in, '*4.csv'))
dirlist = file1csv + file2csv + file3csv + file4csv
for file_csv in dirlist:
  os.remove(file_csv)

print '#### Done removing all unecessary files (*1.csv, *2.csv, *3.csv and *4.csv) ####'

# move the file to a specified directory TODO
#for file_csv in dirlist
#os.rename(os.path.join(path_in), os.path.join(path_out,))

del file1csv, file2csv, file3csv, file4csv, dirlist