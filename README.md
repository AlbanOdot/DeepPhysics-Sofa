# DeepPhysicsSofa

## Ubuntu 
- Eigen lib :       sudo apt install libeigen3-dev
- Python 3 :        sudo apt install python3-dev
- pybind11 :        sudo apt install pybind11-dev
- OpenMP :          sudo apt install libomp-dev
- Intel MKL :       sudo apt install libmkl-full-dev
- SOFA Framework :  https://github.com/sofa-framework/sofa
- SofaPython3 :     https://github.com/sofa-framework/SofaPython3
- Caribou :         https://github.com/jnbrunet

Install the libraries in the given order.

Once project is compiled and installed run :
export DeepPhysicsSofa_ROOT="${PWD}/install"   from build folder
ln -sFfv $(find $DeepPhysicsSofa_ROOT/lib/python3/site-packages -maxdepth 1 -mindepth 1) $(python3 -m site --user-site)
