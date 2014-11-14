### BP Make File
# This makes program "doBP_FULL" which performs the entire pipeline
# from data input to Decimation and through gradient descent.
# The best tutorial so far: http://kiwi.atmos.colostate.edu/fortran/docs/class10.pdf


### NOTES
# doBPhg depends on bp_FULL *.f
# first compile all modules
# then compile the main program
# then move relevant .o files to UnitTesting/
# finally, delete unused .o files

# FC is the symbol for modern fortran macros, including gfortran
# Output flags for FC are:
#    * FC
#	 * FCFLAGS
#	 * FCLIBS

# The first target (usually the executable) is the one that will be created when running the MakeFile 
# To customize the make command, use "make -f <myMakeFileName>"
# Variable names in makefiles are case-sensitive



##### Variables
FC = gfortran
AR = ar
FCFLAGS = -g
FCFLAGS += -I/usr/include		#we have no include files, but this line is good to have in case we want to add them
FCFLAGS = -O14					#Use option 14 for gfortran compiler 
PROGRAMS = GradientDescent doBP_FULL			#full pipeline of 

CO = rm

# Declare gd module targets and .o objects
gdobjects_module = GradientDescent\modarguments.o GradientDescent\modexp.o GradientDescent\modmodel.o GradientDescent\modnonzero.o GradientDescent\modode.o GradientDescent\modoptions.o 

gdobjects_all = $(gdobjects_module) GradientDescent\dlsode.o GradientDescent\exp.o GradientDescent\fit_static_2.o GradientDescent\main.o GradientDescent\model.o GradientDescent\prep_cg.o GradientDescent\prep_data.o GradientDescent\random.o GradientDescent\sim.o

# Declare the object files
objects = bp_mod_CONSTANTS.o bp_mod_BP.o bp_mod_DECIMATION.o bp_allocation.o bp_BeliefPropagation.o bp_calculations.o bp_Decimation.o bp_read_write.o bp_main.o testInd2Sub.o

# Make----------------------------------------------------------------------------------
all: $(PROGRAMS)	#make the executables listed in the 'PROGRAMS' variable

GradientDescent: $(gdobjects_all)
	$(AR) r GradientDescent\libPineda.a $(gdobjects_all)

GradientDescent\modarguments.o:GradientDescent\modarguments.f
	cd GradientDescent; $(FC) -c modarguments.f

GradientDescent\modexp.o:GradientDescent\modexp.f
	cd GradientDescent; $(FC) -c modexp.f

GradientDescent\modmodel.o:GradientDescent\modmodel.f
	cd GradientDescent; $(FC) -c modmodel.f

GradientDescent\modnonzero.o:GradientDescent\modnonzero.f
	cd GradientDescent; $(FC) -c modnonzero.f
	
GradientDescent\modode.o:GradientDescent\modode.f
	cd GradientDescent; $(FC) -c modode.f

GradientDescent\modoptions.o:GradientDescent\modoptions.f
	cd GradientDescent; $(FC) -c modoptions.f

GradientDescent\dlsode.o:GradientDescent\dlsode.f
	cd GradientDescent; $(FC) -c dlsode.f

GradientDescent\exp.o:GradientDescent\exp.f $(gdobjects_module)
	cd GradientDescent; $(FC) -c exp.f

GradientDescent\fit_static_2.o:GradientDescent\fit_static_2.f $(gdobjects_module)
	cd GradientDescent; $(FC) -c fit_static_2.f

GradientDescent\main.o:GradientDescent\main.f $(gdobjects_module)
	cd GradientDescent; $(FC) -c main.f

GradientDescent\model.o:GradientDescent\model.f $(gdobjects_module)
	cd GradientDescent; $(FC) -c model.f

GradientDescen\prep_cg.o:GradientDescent\prep_cg.f $(gdobjects_module)
	cd GradientDescent; $(FC)  -c prep_cg.f

GradientDescent\prep_data.o:GradientDescent\prep_data.f $(gdobjects_module)
	cd GradientDescent; $(FC)  -c prep_data.f

GradientDescent\random.o:GradientDescent\random.f $(gdobjects_module)
	cd GradientDescent; $(FC)  -c random.f

GradientDescent\sim.o:GradientDescent\sim.f $(gdobjects_module)
	cd GradientDescent; $(FC)  -c sim.f

doBP_FULL: $(objects)	#command of how to make 'doBP_FULL' (listed in 'PROGRAMS')
	$(FC) -o doBP_FULL $(FCFLAGS) $(objects) -L.\GradientDescent -lPineda

%bp%.o:%bp%.f
	$(FC) $(FCFLAGS) -c $(objects)

	
clean_dev: 
	mv bp_allocation.o bp_BeliefPropagation.o bp_calculations.o bp_Decimation.o bp_read_write.o *.mod bp_mod_BP.o bp_mod_CONSTANTS.o bp_mod_DECIMATION.o testInd2Sub.o UnitTesting\ ; rm bp_main.o 

clean_usr:
	rm *.o ; rm *.mod
	
#----------------------------------------------------------------------------------------
