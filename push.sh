#!/bin/bash

odir=/media/jerry/internal/data/release/MesoscopicInteractions
c=0

for file in `find .`
do
    if [ $c != 0 ]
    then
        #cmd [option] "$file" >> results.out
        cdir=`pwd`
        echo rsync -avP --exclude . --exclude _Export-example_1.txt --exclude downloads --exclude downloads.tar.gz $file $odir/$file
        rsync -avP --exclude . --exclude _Export-example_1.txt --exclude downloads --exclude downloads.tar.gz $file $odir/$file
        echo cd $odir
        cd $odir
        echo git add *
        git add *
        echo git commit -m "9SEP2020"
        git commit -m "9SEP2020"
        echo git push origin master
        git push origin master
        echo cd $cdir
        cd $cdir
    fi
    c=$((c+1))
done
