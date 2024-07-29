# Certifiable solver for real-time N-view triangulation

This repository contains the code for the fast N-view triangulation solver explained in [the following paper](https://ieeexplore.ieee.org/document/10044919) where the code is publicly available in its [official GitHub page](https://github.com/mergarsal/FastNViewTriangulation).
```BibTeX
@article{garcia2023certifiable,
  title={Certifiable solver for real-time N-view triangulation},
  author={Garcia-Salguero, Mercedes and Gonzalez-Jimenez, Javier},
  journal={IEEE Robotics and Automation Letters},
  volume={8},
  number={4},
  pages={1999--2005},
  year={2023},
  publisher={IEEE}
}
```
## Contributors:
[Mercedes Garcia-Salguero](https://mapir.isa.uma.es/mapirwebsite/?p=1718) <br />
[Javier Gonzalez-Jimenez](https://mapir.isa.uma.es/mapirwebsite/?p=1536) <br />
[Chiang-Heng Chien](https://github.com/C-H-Chien)

## Dependencies
- Eigen (mandatory; tested using version 3.2 and above)
- MATLAB (optional; used for organizing input data for the fast N-view triangulation solver)
- Ceres (optional; used in the example and test. Remove the peice of ceres-solver code if necessary)

## Build for running C++ only
Following the standard build and compile process:
```bash
mkdir build & cd build 
cmake .. 
make -j
```
The compiled examples should be inside the `bin` directory. An example code for generating synthetic camera poses and point correspondences can be found in
```bash
./bin/example_base
```

## Build for running C++ as a solver for MATLAB interface
A matlab example code is provided to give synthetic camera poses and feature point correspondences in which they are fed to the fast N-view triangultion sovler in C++. <br />
Compile the ```mexFunction``` in MATLAB's command window:
```bash
mex examples/fast_multiview_triangulation_mex.cpp src/NViewsCertifier.cpp src/NViewsClass.cpp src/NViewsUtils.cpp utils/generatePointCloud.cpp -I/usr/include/eigen3/
```
Make sure that all the files are under MATLAB's paths. The directory ```-I/usr/include/eigen3/``` links the eigen library, enabling MATLAB to compile the C++ code and directly use ```fast_multiview_triangulation_mex``` function to run N-view triangulation. Refer to ```test_Fast_Nview_Triangulation.m``` for more information.


## Install 
In `build` folder: 
```bash
make install
```
We also provide the uninstall script: 
```bash
make uninstall
```
## Use in external project 
1. Install our library with 
```bash
sudo make install 
```

2. In your project, find the library. 
```bash
find_package(NViewsTrian REQUIRED)
```

3. Include the library as target, e.g., 
```bash
add_executable(example_base ${CMAKE_CURRENT_SOURCE_DIR}/example_base.cpp)
target_link_libraries(example_base NViewsTrian)
```              
