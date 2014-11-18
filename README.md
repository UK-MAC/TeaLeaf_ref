# TeaLeaf

## Compling

In many cases just typing `make` in the required software version will work. This is the case when the mpif90 and mpicc scripts point to the correct compilers.

If the MPI compilers have different names then the build process needs to 
notified of this by defining two environment variables, `MPI_COMPILER` and 
`C_MPI_COMPILER`. 

For example on some Intel systems:

`make MPI_COMPILER=mpiifort C_MPI_COMPILER=mpiicc`

Or on Cray systems:

`make MPI_COMPILER=ftn C_MPI_COMPILER=cc`

Or on GNU:

`make MPI_COMPILER=gfortran C_MPI_COMPILER=gcc`

### OpenMP Build

All compilers use different arguments to invoke OpenMP compilation. A simple 
call to make will invoke the compiler with -O3. This does not usually include 
OpenMP by default. To build for OpenMP for a specific compiler a further 
variable must be defined, `COMPILER` that will then select the correct option 
for OpenMP compilation. 

For example with the Intel compiler:

`make COMPILER=INTEL`

Which then append the -openmp to the build flags.

Other supported compiler that will be recognise are:-

* CRAY
* SUN
* GNU
* IBM
* PATHSCALE
* PGI

The default flags for each of these is show below:-

* INTEL: -O3 -ipo
* SUN: -fast
* GNU: -ipo
* XL: -O5
* PATHSCLE: -O3
* PGI: -O3 -Minline
* CRAY: -em  Note: that by default the Cray compiler with pick the optimum 
options for performance.

### Other Flags

The default compilation with the COMPILER flag set chooses the optimal 
performing set of flags for the specified compiler, but with no hardware 
specific options or IEEE compatability.

To produce a version that has IEEE compatiblity a further flag has to be set on 
the compiler line.

`make COMPILER=INTEL IEEE=1`

This flag has no effect if the compiler flag is not set because IEEE options 
are always compiler specific.

For each compiler the flags associated with IEEE are shown below:-

* INTEL: -fp-model strict –fp-model source –prec-div –prec-sqrt
* CRAY: -hpflex_mp=intolerant
* SUN: -fsimple=0 –fns=no
* GNU: -ffloat-store
* PGI: -Kieee
* PATHSCALE: -mieee-fp
* XL: -qstrict –qfloat=nomaf

Note that the MPI communications have been written to ensure bitwise identical 
answers independent of core count. However under some compilers this is not 
true unless the IEEE flags is set to be true. This is certainly true of the 
Intel and Cray compiler. Even with the IEEE options set, this is not guarantee 
that different compilers or platforms will produce the same answers. Indeed a 
Fortran run can give different answers from a C run with the same compiler, 
same options and same hardware.

Extra options can be added without modifying the makefile by adding two further 
flags, `OPTIONS` and `C_OPTIONS`, one for the Fortran and one for the C options.

`make COMPILER=INTEL OPTIONS=-xavx C_OPTIONS=-xavx`

Finally, a `DEBUG` flag can be set to use debug options for a specific compiler.

`make COMPILER=PGI DEBUG=1`

These flags are also compiler specific, and so will depend on the `COMPILER` 
environment variable.

So on a system without the standard MPI wrappers, for a build that requires 
OpenMP, IEEE and AVX this would look like so:-


make COMPILER=INTEL MPI_COMPILER=mpiifort C_MPI_COMPILER=mpiicc IEEE=1 \
OPTIONS="-xavx" C_OPTIONS="-xavx"

### File Input

The contents of tea.in defines the geometric and run time information, apart from task and thread counts.

A complete list of options is given below, where `<R>` shows the option takes a real number as an argument. Similarly `<I>` is an integer argument.

`initial_timestep <R>`

Set the initial time step for TeaLeaf. This time step stays constant through the entire simulation. The default value is 

`end_time <R>`

Sets the end time for the simulation. When the simulation time is greater than this number the simulation will stop.

`end_step <I>`

Sets the end step for the simulation. When the simulation step is equal to this then simulation will stop.

In the event that both the above options are set, the simulation will terminate on whichever completes first.

`xmin <R>`

`xmax <R>`

`ymin <R>`

`ymax <R>`

The above four options set the size of the computational domain. The default domain size is a 10cm square. 

`x_cells <I>`

`y_cells <I>`

The two options above set the cell count for each coordinate direction. The default is 10 cells in each direction.

The geometric information and initial conditions are set using the following keywords with three possible variations. Note that state 1 is always the ambient material and any geometry information is ignored. Areas not covered by other defined states receive the energy and density of state 1.

`state <I> density <R> energy <R> geometry rectangle xmin <R> ymin <R> xmax <R> ymax <R> `

Defines a rectangular region of the domain with the specified energy and density.

`state <I> density <R> energy <R> geometry circle xmin <R> ymin <R> radius <R>`

Defines a circular region of the domain with the specified energy and density.

`state <I> density <R> energy <R> geometry point xmin <R> ymin <R>`

Defines a cell in the domain with the specified energy and density.

Note that the generator is simple and the defined state completely fills a cell with which it intersects. In the case of over lapping regions, the last state takes priority. Hence a circular region will have a stepped interface and a point data will fill the cell it lies in with its defined energy and density.

`visit_frequency <I>`

This is the step frequency of visualisations dumps. The files produced are text base VTK files and are easily viewed in an application such as ViSit. The default is to output no graphical data. Note that the overhead of output is high, so should not be invoked when performance benchmarking is being carried out.

`summary_frequency <I>`

This is the step frequency of summary dumps. This requires a global reduction and associated synchronisation, so performance will be slightly affected as the frequency is increased. The default is for a summary dump to be produced every 10 steps and at the end of the simulation.

`tl_ch_cg_presteps  <I>`

This option specifies the number of Conjugate Gradient iterations completed before the Chebyshev method is started. This is necessary to provide approximate minimum and maximum eigen values to start the Chebyshev method. The default value is 30.

`tl_ppcg_inner_steps <I>`

The default value is 10.

`tl_ch_cg_epslim`

The default value is 1e-5.

`tl_check_result`

The default for this option is off.

`tl_ch_cg_errswitch`

The default for this is off.

`tl_preconditioner_on`

This keyword invokes the pre-conditioner. The only pre-conditioner available is a diagonal scaling. This is a simple precoditioner and may not accelerate the time to solution or reduce the iteration count. By default, no pre-conditioner is used.

`tl_use_jacobi`

This keyword selects the Jacobi method to solve the linear system. Note that this a very slowly converging method compared to other options. This is the default method is no method is explicitly selected.

`tl_use_cg`

This keyword selects the Conjugate Gradient method to solve the linear system.

`tl_use_ppcg`

This keyword selects the Conjugate Gradient method to solve the linear system.

`tl_use_chebyshev`

This keyword selects the Chebyshev method to solve the linear system.

`profiler_on`

This option turns the code's coarse grained internal profiler end. Timing information is reported at the end of the simulation in the tea.out file. The default is no profiling.

`tl_max_iters <I>`

This option provides an upper limit of the number of iterations used for the linear solve in a step. If this limit is reached, then the solution vector at this iteration is used as the solution, even if the convergence criteria has not been met. For this reason, care should be taken in the comparison of the performance of a slowly converging method, such as Jacobi, as the convergence criteria may not have been met for some of the steps. The default value is 1000.

`tl_eps <R>`

This option sets the convergence criteria for the selected solver. It uses a least squares measure of the residual. The default value is 1.0e-10.

`tl_coefficient_density

This option uses the density as the conduction coefficient. This is the default option.

`tl_coefficient_inverrse_density

This option uses the inverse density as the conduction coefficient.

`test_problem <I>`

This keyword selects a standard test with a "known" solution. Test problem 1 is automatically generated if the tea.in file does not exist. Test problems 2-5 are shipped in the TeaLeaf repository. Note that the known solution for an iterative solver is not an analytic solution but is the solution for a single core simulation with IEEE options enabled with the Intel compiler and a strict convergence of 1.0e-15. The difference to the expected solution is reported at the end of the simulation in the tea.out file. There is no default value for this option.



