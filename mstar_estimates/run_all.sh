#!/bin/bash

python mstar.SEDs.general.py vandenberk Hopkins07 Full
python mstar.SEDs.general.py vandenberk Shen20 A
python mstar.SEDs.general.py vandenberk Shen20 B

python mstar.SEDs.general.py Richards06 Hopkins07 Full
python mstar.SEDs.general.py Richards06 Shen20 A
python mstar.SEDs.general.py Richards06 Shen20 B

python mstar.lrt.general.py Hopkins07 Full
python mstar.lrt.general.py Shen20 A
python mstar.lrt.general.py Shen20 B
