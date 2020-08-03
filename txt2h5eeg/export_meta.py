import os

def main():
    print('--- Export Meta ---')

    # check meta folder
    if (not os.path.isdir('./meta')):
        os.system('mkdir meta')
        print('[*] ./meta created.')
    else:
        print('[*] ./meta exists.')

    # check isaux.txt
    if (not os.path.isfile('./meta/isaux.txt')):
        print('[*] ./meta/isaux.txt missing.')
        cond_isaux = False
    else:
        cond_isaux = True
        print('[*] ./meta/isaux.txt exists.')

    # check aux_labels.txt
    if (not os.path.isfile('./meta/aux_labels.txt')):
        print('[*] ./meta/aux_labels.txt missing.')
        cond_aux_labels = False
    else:
        cond_aux_labels = True
        print('[*] ./meta/aux_labels.txt exists.')

    # check channel_labels.txt
    if (not os.path.isfile('./meta/channel_labels.txt')):
        print('[*] ./meta/channel_labels.txt missing.')
        cond_channel_labels = False
    else:
        cond_channel_labels = True
        print('[*] ./meta/channel_labels.txt exists.')

    # export meta
    cond_exportmeta = (not cond_isaux) or (not cond_aux_labels) or (not cond_channel_labels)
    if (cond_exportmeta):
        print('[*] Exporting metadata..')

        # Get list of existing .txt files
        Fnames = []
        for fi in os.listdir('./'):
            if fi.endswith('.txt'):
                Fnames.append(fi)
                print('\t' + fi)
        if (len(Fnames) > 1):
            print('[*] Using first .txt file for meta export.')

        # Read Neuroworks .txt
        c = 0
        with open(Fnames[0]) as infile:
            for l in infile:
                # Channel labels line
                if ('C001' in l):
                    # Parse labels
                    Clabels = []
                    for item in (l[1:].split()):
                        if item.startswith('C'):
                            Clabels.append(item)
                    #print(Clabels)
                    print('(*) Found ' + str(len(Clabels)) + ' channel labels.')

                    # Write
                    ofile_isaux = open('./meta/isaux.txt','w')
                    ofile_aux_labels = open('./meta/aux_labels.txt','w')
                    ofile_channel_labels = open('./meta/channel_labels.txt','w')
                    for i in range(len(Clabels)):
                        print('0',file=ofile_isaux)
                        print(Clabels[i],file=ofile_channel_labels)
                    ofile_isaux.close()
                    ofile_aux_labels.close()
                    ofile_channel_labels.close()
                if (c > 18):
                    #print(l[1:].split())
                    #print(len(l[1:].split()))
                    break
                c = c + 1
        print('[*] Done.')
    else:
        print('[*] Existing metadata found, exiting.')

main()
