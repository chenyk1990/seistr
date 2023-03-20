**seistr**
======

## Description

**seistr** is an open-source Matlab package for 2D/3D local slope estimation and structural filtering. The seistr package has a variety of applications in both exploration and earthquake seismology, including but not limited to seismic denoising, seismic reconstruction, seismic diffraction separation, constrained LSRTM, constrained FWI, etc. To the best of our knowledge, this is the first (and only) open-source Matlab package for calculating the seismic slope of an input seismic gather/cube/image and performing structure-oriented filtering. 

Since March, 19, 2023, this official site of Matlab-version seistr package has been moved to https://github.com/aaspip/MATseistr. The current site of the package is no longer maintained.

The python version of the seistr package can be found at https://github.com/aaspip/pyseistr.

## Reference
    Wang, H., Chen, Y., Saad, O.M., Chen, W., Oboué, Y.A.S.I., Yang, L., Fomel, S. and Chen, Y., 2022. A Matlab code package for 2D/3D local slope estimation and structural filtering. Geophysics, 87(3), pp.F1–F14.
    
BibTeX:

	@article{seistr,
	  title={A Matlab code package for 2{D}/3{D} local slope estimation and structural filtering},
	  author={Hang Wang and Yunfeng Chen and Omar M. Saad and Yapo Abol\'{e} Serge Innocent Obou\'{e} and Liuqing Yang and Sergey Fomel and Yangkang Chen},
	  journal={Geophysics},
	  volume={87},
	  number={3},
	  issue={3},
	  pages={F1–F14},
	  year={2022},
	  doi={10.1190/geo2021.0266.1},
	  publisher={Society of Exploration Geophysicists}
	}

-----------
## Copyright
    seistr developing team, 2021-present
-----------

## License
    GNU General Public License, Version 3
    (http://www.gnu.org/copyleft/gpl.html)   

-----------

## Install
Using the latest version

    git clone https://github.com/chenyk1990/seistr
    cd seistr
    addpath(genpath('./')); #in Matlab command line
    
-----------
## Examples
    The "demo" directory contains all runable scripts to demonstrate different applications of seistr. 

-----------
## Dependence Packages
* Matlab 2015 and later versions
    
-----------
## Development
    The development team welcomes voluntary contributions from any open-source enthusiast. 
    If you want to make contribution to this project, feel free to contact the development team. 

-----------
## Contact
    Regarding any questions, bugs, developments, collaborations, please contact  
    Yangkang Chen
    chenyk2016@gmail.com

-----------
## Gallery
The gallery figures of the seistr package can be found at
    https://github.com/chenyk1990/gallery/tree/main/seistr
Each figure in the gallery directory corresponds to a DEMO script in the "demo" directory with the exactly the same file name.

The first example is a 2D synthetic test generated by test_seistr_syn2d.m in the demos directory.
<img src='https://github.com/chenyk1990/gallery/blob/main/seistr/test_seistr_syn2d.png' alt='Syn2D' width=960/>

The second example is a 3D real data test generated by test_seistr_real3d.m in the demos directory.
<img src='https://github.com/chenyk1990/gallery/blob/main/seistr/test_seistr_real3d.png' alt='Real3D' width=960/>

