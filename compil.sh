#gfortran -C -traceback -o interp3d bspline_module.f90 bspline_sub_module.f90 bspline_oo_module.f90 test_interp3d.f90
#gfortran  -o interp3d test_interp3d.f90 bspline_module.o bspline_sub_module.o bspline_oo_module.o


gfortran -ffpe-summary='none' -w -ffree-line-length-none -c interp.f90
gfortran -ffpe-summary='none' -w -ffree-line-length-none -c bspline_sub_module.f90 
gfortran -ffpe-summary='none' -w -ffree-line-length-none -c bspline_oo_module.f90 
gfortran -ffpe-summary='none' -w -ffree-line-length-none -c bspline_module.f90 
gfortran -ffpe-summary='none' -w -ffree-line-length-none -c distrib.f90 
gfortran -ffpe-summary='none' -w -ffree-line-length-none -c misc.f90 
gfortran -ffpe-summary='none' -w -ffree-line-length-none -o dyn dynamic.f90 dynlib.f90 bspline_module.o bspline_sub_module.o bspline_oo_module.o distrib.o misc.o interp.o

#ifort -c interp.f90
#ifort -c bspline_sub_module.f90 
#ifort -c bspline_oo_module.f90 
#ifort -c bspline_module.f90 
#ifort -c distrib.f90 
#ifort -c misc.f90 
#ifort -o dyn dynamic.f90 dynlib.f90 bspline_module.o bspline_sub_module.o bspline_oo_module.o distrib.o misc.o interp.o
