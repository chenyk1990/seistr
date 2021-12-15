This is the open-source Matlab package for 2D/3D local slope estimation and structural filtering.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
In the "code2D" folder, there is a 2D example. Its main script is "test_2D_SOF.m" which calls all subfunctions in this folder. The script has three parts including:

1. generating the test data;
2. 2D local slope estimation;
3. the 2D structural filtering.

To get the 2D demo, please run the "test_2D_SOF.m". 

The user can replace the data generation part by loading existing synthetic or field data file to test their own dataset.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
In the "code3D" folder, there is a 3D demo. Its main script is "test_3D_SOF.m" which also calls all subfunctions in this folder. The script also has three parts including:

1. loading the 3D data "real3d.bin";
2. 3D local slope estimation;
3. the 3D structural filtering.

To get the 3D demo, please run the "test_3D_SOF.m". 

The user can obtain the input 3D source data file "real3d.bin" at https://github.com/chenyk1990/reproducible_research/blob/master/drr3d/matfun/real3d.bin.

The user can also replace the data loading part with their own 3D dataset.



