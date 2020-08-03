# Laste edited on January 7, 2017
# Jiarui Wang :: jwang04@g.harvard.edu

import os

# Color class
class col:
    H = '\033[95m' # Header
    B = '\033[94m' # Blue
    G = '\033[92m' # Green
    W = '\033[93m' # Warning
    E = '\033[91m' # Error
    e = '\033[0m' # ENDC

# Converter class
class converter:

    #
    # Constructor
    #
    # Stores exported and finished directories
    def __init__(self, Dirs, pid):
        self.exDir = Dirs[0]
        self.fiDir = Dirs[1]
        self.pid = pid
        self.exportList = []

    #---------------------------------------------------------------------------
    # "HIDDEN" FUNCTIONS
    #
    # Don't call these directly
    #---------------------------------------------------------------------------
    
    #
    # Check that there are files to convert: self.exportList is empty
    #
    def checkConvert(self):
        if (len(self.exportList) == 0):
            print(col.E+'[!] Error: no files to convert.'+col.e)
            exit()

    # 
    # Makes sure there are output directories in self.fiDir
    # 
    def involveOutput(self):
        os.system('mkdir '+self.fiDir+'/'+self.pid)
        os.system('mkdir '+self.fiDir+'/'+self.pid+'/tmp')

    #
    # Prints and executes commands on the shell
    #
    def shell(self, command):
        print(command)
        os.system(command)

    #
    # Checks whether or not given filename is the matching annotation file
    #
    def isAnnot(self, efn, afn):
        afn = afn.split('.ent.txt')[0]
        afn = afn.split('-')[0]
        efn = efn.split('.txt')[0]
        efn = efn.split('-')[1]
        #print(afn)
        #print(efn)
        if (afn == efn):
            return True
        elif (afn.upper() == efn.upper()):
            return True
        else:
            return False

    #---------------------------------------------------------------------------
    # PUBLIC FUNCTIONS
    #---------------------------------------------------------------------------

    #
    # Make list of exported text file names under the exported directory
    #
    # fills self.exportList
    #
    def buildExportList(self):
        # Make self.exportList (initialized in constructor)
        dirList = next(os.walk(self.exDir))[1]
        for d in dirList:
            txtList = next(os.walk(self.exDir+'/'+d))[2]
            # Only add exported text files for given patient ID
            if (self.pid == d):
                for ti in txtList:
                    #self.exportList.append(self.exDir+'/'+self.pid+'/'+ti)
                    self.exportList.append(ti)
        #print(self.exportList)
        
        # Avoid .DS_Store left by macOS finder
        if ('.DS_Store' in self.exportList):
            self.exportList.remove('.DS_Store')

        # Check to make sure txt files were found
        if (len(self.exportList) == 0):
            print(col.E+'[!] Error: no files found for patient ID'+col.e)
            exit()
        # If files exist, display them
        else:
            print('[!] Patient '+self.pid+' has the following files:')
            for ename in self.exportList:
                print('\t'+ename)
                #os.system('tail '+ename) # Make sure there's actually data
   
    #
    # Entry point for main conversion
    #
    def chunkExported(self):
        # Make output folders if they don't already exist
        self.involveOutput()
        # Make sure buildExportList has been called
        self.checkConvert()

        # Define locations
        tmpDir = self.fiDir + '/' + self.pid + '/tmp'
        ffDir = self.fiDir + '/' + self.pid
        efDir = self.exDir + '/' + self.pid
        currDir = os.getcwd()
        #pDir = self.fiDir + '/' + self.pid

        # self.exportList main loop
        for efile in self.exportList:
            # Show conversion start message
            print('\n[!] Beginning conversion for: '+efile)
            # Bash friendly filename
            Efile = efile.replace(' ','\ ')
            # Move text file to tmp directory
            self.shell('mv '+efDir+'/'+Efile+' '+tmpDir+'/')
            # Clean intermediate files
            self.shell('rm '+tmpDir+'/AA_00000000_*')
            #Make symbolic link in tmp directory
            #self.shell('ln -s '+efDir+'/'+Efile+' '+tmpDir+'/'+Efile)

            # --- Chunking ---
            self.shell('cp '+currDir+'/'+'chunkraw.py'+' '+tmpDir+'/')
            os.chdir(tmpDir)
            self.shell('python3 chunkraw.py '+Efile)
            os.chdir(currDir)
            self.shell('rm '+tmpDir+'/chunkraw.py')

            # Move text file back to export directory
            self.shell('mv '+tmpDir+'/'+Efile+' '+efDir+'/')

            # --- Removing OFF ---
            self.shell('cp '+currDir+'/'+'removeOff.py'+' '+tmpDir+'/')
            os.chdir(tmpDir)
            self.shell('python3 removeOff.py')
            os.chdir(currDir)
            self.shell('rm '+tmpDir+'/removeOff.py')
            
            # --- Peter ---
            self.shell('cp '+currDir+'/'+'peter.py'+' '+tmpDir+'/')
            os.chdir(tmpDir)
            self.shell('python3 peter.py '+currDir+'/'+efDir+'/'+Efile)
            os.chdir(currDir)
            self.shell('rm '+tmpDir+'/peter.py')

            # --- Mainfilter ---
            self.shell('cp '+currDir+'/'+'*.m'+' '+tmpDir+'/')
            os.chdir(tmpDir)
            self.shell('matlab -nodesktop -nosplash -r mainfilter')
            # save preprocessing outputs
            self.shell('mv preprocessing_output PLI_'+Efile.replace('.txt',''))
            os.chdir(currDir)
            self.shell('rm '+tmpDir+'/*.m')

            # --- Delete .csv files to save space ---
            #self.shell('rm '+tmpDir+'/*.csv')

            # --- Export meta ---
            self.shell('cp '+currDir+'/'+'export_meta.py'+' '+efDir+'/')
            os.chdir(efDir)
            self.shell('python3 export_meta.py')
            os.chdir(currDir)
            self.shell('rm '+efDir+'/export_meta.py')

            # --- Convert2H5eeg ---
            # Find annotation files
            #os.chdir(efDir+'/meta')
            metaList = next(os.walk(efDir+'/meta'))[2]
            afile = ''
            for afilen in metaList:
                if (self.isAnnot(efile,afilen)):
                    afile = afilen
                    break

            metaErrors = 0
            # Find annotation file
            if (afile == ''):
                print(col.W+'[!] Warning: annotation file for '+efile+\
                    ' not found. Using null template'+col.e)
                afile = 'annotations.txt.null'
                self.shell('cp '+currDir+'/'+afile+' '+efDir+'/meta/')
                #metaErrors = metaErrors + 1
            else:
                print('[!] Annotations: using file '+afile)

            # Check for channel_labels.txt
            if ('channel_labels.txt' in metaList):
                print('[!] channel_labels.txt: found!')
            else:
                print(col.E+'[!] Error: channel_labels.txt not found.'+col.e)
                metaErrors = metaErrors + 1

            # Check for aux_labels.txt
            if ('aux_labels.txt' in metaList):
                print('[!] aux_labels.txt: found!')
            else:
                print(col.E+'[!] Error: aux_labels.txt not found.'+col.e)
                metaErrors = metaErrors + 1

            # Check for isaux.txt
            if ('isaux.txt' in metaList):
                print('[!] isaux.txt: found!')
            else:
                print(col.E+'[!] Error: isaux.txt not found.'+col.e)
                metaErrors = metaErrors + 1

            # Print instructions for any missing files 
            if (metaErrors > 0):
                me = str(metaErrors)
                print(col.H+'[!] '+me+' metafiles not found. Add these'+\
                    ' files to meta/ within [exported]/[patientID]/'+col.e)
                exit()

            # Copy meta files to finished directory
            Afile = afile.replace(' ','\ ')
            self.shell('cp '+efDir+'/meta/'+Afile+' '+ffDir+'/annotations.txt')
            self.shell('cp '+efDir+'/meta/isaux.txt'+' '+ffDir+'/')
            self.shell('cp '+efDir+'/meta/aux_labels.txt'+' '+ffDir+'/')
            self.shell('cp '+efDir+'/meta/channel_labels.txt'+' '+ffDir+'/')
            self.shell('cp '+currDir+'/'+'*.m'+' '+ffDir+'/')
            
            # Execute
            os.chdir(ffDir)
            self.shell('matlab -nodesktop -nosplash -r convert2h5eeg3')
            EfileN = Efile.split('.txt')[0]
            EfileN = EfileN.split('-')[1]
            self.shell('mv '+self.pid+'.hdf5 '+EfileN+'.hdf5')
            os.chdir(currDir)
            
            # Clean up 
            self.shell('rm '+ffDir+'/*.m')
            self.shell('rm '+ffDir+'/channel_labels.txt')
            self.shell('rm '+ffDir+'/aux_labels.txt')
            self.shell('rm '+ffDir+'/isaux.txt')
            self.shell('rm '+ffDir+'/annotations.txt')

            # Clean tmp folder
            self.shell('rm -f '+tmpDir+'/*')
            
            # Show conversion end message
            print(col.G+'[!] Finished conversion for: '+efile+col.e)

    #
    # For debugging
    #
    def printContents(self):
        print(self.exDir)
        print(self.fiDir)

