# DeepPhysicsSofa

## Ubuntu 
- Eigen lib :
```bash
sudo apt install libeigen3-dev
```
- Python 3 :
```bash
sudo apt install python3-dev
```
- pybind11 :
```bash
sudo apt install pybind11-dev
```
- OpenMP :
```bash
sudo apt install libomp-dev
```
- Intel MKL :
```bash
sudo apt install libmkl-full-dev
```
- SOFA Framework :  https://github.com/sofa-framework/sofa
- SofaPython3 :     https://github.com/sofa-framework/SofaPython3
- Caribou :         https://github.com/jnbrunet


Install the libraries in the given order.

Once project is compiled and installed run from build folder :
```bash
export DeepPhysicsSofa_ROOT="$(pwd)/install"
ln -sFfv $(find $DeepPhysicsSofa_ROOT/lib/python3/site-packages -maxdepth 1 -mindepth 1) $(python3 -m site --user-site)
```
