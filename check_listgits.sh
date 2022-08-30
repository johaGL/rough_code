#!/bin/bash
cd $HOME
for i in $(cat $HOME/rough_code/listgits) ; do
    echo "=>     ${i}       <="; 
    cd $i; 
    git status;
    git log | head -n 4;
    cd $HOME;
    echo "-  -  -  -  -  -  -  -  -  --  -  -  -  -  -  -  -  -  -";
    echo "";
done

