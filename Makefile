#Crown Copyright 2014 AWE.
#
# This file is part of TeaLeaf.
#
# TeaLeaf is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the 
# Free Software Foundation, either version 3 of the License, or (at your option) 
# any later version.
#
# TeaLeaf is distributed in the hope that it will be useful, but 
# WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
# FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more 
# details.
#
# You should have received a copy of the GNU General Public License along with 
# TeaLeaf. If not, see http://www.gnu.org/licenses/.

#  @brief Makefile for CloverLeaf
#  @author David Beckingsale, Wayne Gaudin
#  @details Agnostic, platform independent makefile for the TeaLeaf benchmark code.

# Agnostic, platform independent makefile for the TeaLeaf benchmark code.
# It is not meant to be clever in anyway, just a simple build out of the box script.
# Just make sure mpif90 is in your path. It uses mpif90 even for all builds because this abstracts the base
#  name of the compiler. If you are on a system that doesn't use mpif90, just replace mpif90 with the compiler name
#  of choice. The only mpi dependencies in this non-MPI version are mpi_wtime in timer.f90.

# There is no single way of turning OpenMP compilation on with all compilers.
# The known compilers have been added as a variable. By default the make
#  will use no options, which will work on Cray for example, but not on other
#  compilers.
# To select a OpenMP compiler option, do this in the shell before typing make:-
#
#  export COMPILER=INTEL       # to select the Intel flags
#  export COMPILER=SUN         # to select the Sun flags
#  export COMPILER=GNU         # to select the Gnu flags
#  export COMPILER=CRAY        # to select the Cray flags
#  export COMPILER=PGI         # to select the PGI flags
#  export COMPILER=PATHSCALE   # to select the Pathscale flags
#  export COMPILER=XLF         # to select the IBM Xlf flags

# or this works as well:-
#
# make COMPILER=INTEL
# make COMPILER=SUN
# make COMPILER=GNU
# make COMPILER=CRAY
# make COMPILER=PGI
# make COMPILER=PATHSCALE
# make COMPILER=XLF
#

# Don't forget to set the number of threads you want to use, like so
# export OMP_NUM_THREADS=4

# usage: make                     # Will make the binary
#        make clean               # Will clean up the directory
#        make DEBUG=1             # Will select debug options. If a compiler is selected, it will use compiler specific debug options
#        make IEEE=1              # Will select debug options as long as a compiler is selected as well
# e.g. make COMPILER=INTEL MPI_COMPILER=mpiifort C_MPI_COMPILER=mpiicc DEBUG=1 IEEE=1 # will compile with the intel compiler with intel debug and ieee flags included

ifndef COMPILER
  MESSAGE=select a compiler to compile in OpenMP, e.g. make COMPILER=INTEL
endif

OMP_INTEL     = -openmp
OMP_SUN       = -xopenmp=parallel -vpara
OMP_GNU       = -fopenmp
OMP_CRAY      =
OMP_PGI       = -mp=nonuma
OMP_PATHSCALE = -mp
OMP_XLF       = -qsmp=omp -qthreaded
OMP=$(OMP_$(COMPILER))

FLAGS_INTEL     = -O3 -ipo
FLAGS_SUN       = -fast
FLAGS_GNU       = -O3
FLAGS_CRAY      = -em -ra -h acc_model=fast_addr:no_deep_copy:auto_async_all
FLAGS_PGI       = -O3
FLAGS_PATHSCALE = -O3
FLAGS_XLF       = -O3
FLAGS_          = -O3
CFLAGS_INTEL     = -O3 -restrict -fno-alias -ipo
CFLAGS_SUN       = -fast
CFLAGS_GNU       = -O3
CFLAGS_CRAY      = -em -h list=a
CFLAGS_PGI       = -O3
CFLAGS_PATHSCALE = -O3
CFLAGS_XLF       = -O5 -qextname=flush:ideal_gas_kernel_c:viscosity_kernel_c:pdv_kernel_c:revert_kernel_c:accelerate_kernel_c:flux_calc_kernel_c:advec_cell_kernel_c:advec_mom_kernel_c:reset_field_kernel_c
CFLAGS_          = -O3

ifdef DEBUG
  FLAGS_INTEL     = -O0 -g -debug all -check all -traceback -check noarg_temp_created
  FLAGS_SUN       = -O0 -xopenmp=noopt -g
  FLAGS_GNU       = -O0 -g
  FLAGS_CRAY      = -O0 -g -em -eD
  FLAGS_PGI       = -O0 -g -C -Mchkstk -Ktrap=fp -Mchkfpstk
  FLAGS_PATHSCALE = -O0 -g
  FLAGS_XLF       = -O0 -g
  FLAGS_          = -O0 -g
  CFLAGS_INTEL    = -O0 -g -c -debug all -traceback -restrict
  CFLAGS_CRAY     = -O0 -g -em -eD
endif

ifdef IEEE
  I3E_INTEL     = -fp-model strict -fp-model source -prec-div -prec-sqrt
  I3E_SUN       = -fsimple=0 -fns=no
  I3E_GNU       = -ffloat-store
  I3E_CRAY      = -hflex_mp=intolerant
  I3E_PGI       = -Kieee
  I3E_PATHSCALE = -mieee-fp
  I3E_XLF       = -qfloat=nomaf
  I3E=$(I3E_$(COMPILER))
endif

FLAGS=$(FLAGS_$(COMPILER)) $(OMP) $(I3E) $(OPTIONS)
CFLAGS=$(CFLAGS_$(COMPILER)) $(OMP) $(I3E) $(C_OPTIONS) -c
MPI_COMPILER=mpif90
C_MPI_COMPILER=mpicc

tea_leaf: c_lover *.f90 Makefile
	$(MPI_COMPILER) $(FLAGS)	\
	data.f90			\
	definitions.f90			\
	clover.f90			\
	report.f90			\
	timer.f90			\
	parse.f90			\
	read_input.f90			\
	initialise_chunk_kernel.f90	\
	initialise_chunk.f90		\
	build_field.f90			\
	update_halo_kernel.f90		\
	update_halo.f90			\
	ideal_gas_kernel.f90		\
	ideal_gas.f90			\
	start.f90			\
	generate_chunk_kernel.f90	\
	generate_chunk.f90		\
	initialise.f90			\
	field_summary_kernel.f90	\
	field_summary.f90		\
	viscosity_kernel.f90		\
	viscosity.f90			\
	calc_dt_kernel.f90		\
	calc_dt.f90			\
	timestep.f90			\
	accelerate_kernel.f90		\
	accelerate.f90			\
	revert_kernel.f90		\
	revert.f90			\
	PdV_kernel.f90			\
	PdV.f90				\
	flux_calc_kernel.f90		\
	flux_calc.f90			\
	advec_cell_kernel.f90		\
	advec_cell_driver.f90		\
	advec_mom_kernel.f90		\
	advec_mom_driver.f90		\
	advection.f90			\
	reset_field_kernel.f90		\
	reset_field.f90			\
	set_field_kernel.f90            \
	set_field.f90                   \
	tea_leaf_kernel.f90             \
	tea.f90                         \
	hydro.f90		        \
	visit.f90			\
	tea_leaf.f90			\
	accelerate_kernel_c.o           \
	PdV_kernel_c.o                  \
	flux_calc_kernel_c.o            \
	revert_kernel_c.o               \
	reset_field_kernel_c.o          \
	ideal_gas_kernel_c.o            \
	viscosity_kernel_c.o            \
	advec_mom_kernel_c.o            \
	advec_cell_kernel_c.o           \
	tea_leaf_kernel_c.o             \
	-o tea_leaf; echo $(MESSAGE)

c_lover: *.c Makefile
	$(C_MPI_COMPILER) $(CFLAGS)     \
	accelerate_kernel_c.c           \
	PdV_kernel_c.c                  \
	flux_calc_kernel_c.c            \
	revert_kernel_c.c               \
	reset_field_kernel_c.c          \
	ideal_gas_kernel_c.c            \
	viscosity_kernel_c.c            \
	field_summary_kernel_c.c        \
	advec_mom_kernel_c.c            \
	advec_cell_kernel_c.c           \
	tea_leaf_kernel_c.c


clean:
	rm -f *.o *.mod *genmod* tea_leaf tea.in.tmp tea.out
