# CHERAB
sudo apt-get update//
sudo apt install python3-pip
pip install -U cython==3.0a5
pip install numpy
pip install scipy
pip install matplotlip
python3 setup.py install --user in raysect-source/source
python3 setup.py install --user in cherab-source/core
cherab-source/solps should be installed from the official repository
once it is installed do in cherab-source/solps/:
python3 setup.py install --user

pip install mendeleev
