#gfortran -C -traceback -o interp3d bspline_module.f90 bspline_sub_module.f90 bspline_oo_module.f90 test_interp3d.f90
#gfortran  -o interp3d test_interp3d.f90 bspline_module.o bspline_sub_module.o bspline_oo_module.o


gfortran -c interp.f90
gfortran -c bspline_module.f90
gfortran -c bspline_sub_module.f90
gfortran -c bspline_oo_module.f90
gfortran -c distrib.f90
gfortran -c misc.f90
gfortran -o dyn dynamic.f90 dynlib.f90 bspline_module.o bspline_sub_module.o bspline_oo_module.o distrib.o misc.o 
