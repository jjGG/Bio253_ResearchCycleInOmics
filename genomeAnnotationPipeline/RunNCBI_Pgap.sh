#!/bin/sh

#  RunNCBI_Pgap.sh
#  
#
#  Created by jonasgrossmann on 22.11.2024.
#  


# from https://github.com/ncbi/pgap/wiki/Quick-Start
wget -O pgap.py https://github.com/ncbi/pgap/raw/prod/scripts/pgap.py

# make it executable
chmod +x pgap.py

# download docker and all resources necessairy (quite some GB)
./pgap.py --update

# how to invoke my own SA genome?
# ./pgap.py -r -o <results> -g <fasta> -s 'Staphylococcus aureus'
./pgap.py -r -o -g -s 'Staphylococcus aureus'


