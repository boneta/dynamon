CC    = gcc
FC    = gfortran

FLGS  = -I. -J. -w
FLGS += -finit-local-zero -frecord-marker=4
FLGS += -O3

# FLGS += -DHACK_CHRG  ## Code modified to be always activated and not depend on\
                          compilation flags that forces to choose either semiemp or abinitio \
                          What is the f** purpose of SKIP_ABINITIO variable then??????? (╯°□°）╯︵ ┻━┻
FLGS += -DHACK_LAMB
FLGS += -DHACK_KIE

OBJS =	BLASD.o\
		LAPACK_DLAMC1.o\
		LAPACK.o\
		boundary.o\
		definitions.o\
		constants.o\
		elements.o\
		io_units.o\
		printing.o\
		string.o\
		time.o\
		files.o\
		parsing.o\
		conjugate_gradient.o\
		diagonalization.o\
		linear_algebra.o\
		linear_equations.o\
		random.o\
		sort.o\
		special_functions.o\
		statistics.o\
		baker_optimization.o\
		atoms.o\
		sequence.o\
		symmetry.o\
		atom_manipulation.o\
		connectivity.o\
		geometry.o\
		sas.o\
		transformation.o\
		superimpose.o\
		coordinate_io.o\
		pdb_io.o\
		xyz_io.o\
		zmatrix.o\
		zmatrix_io.o\
		constraint.o\
		hack_spline.o\
		wham.o\
		hack_pert.o\
		mm_file_data.o\
		mm_file_io.o\
		mm_terms.o\
		mm_system.o\
		mm_system_edit.o\
		mm_system_io.o\
		build_coordinates.o\
		mopac_data.o\
		mopac_parameters.o\
		energy_covalent.o\
		energy_nb.o\
		hessian_update.o\
		gaussian_basis.o\
		mopac_density.o\
		mopac_fock.o\
		mopac_hamiltonian.o\
		mopac_integrals.o\
		mopac_scf.o\
		mopac_analysis.o\
		mopac_gradients.o\
		quantum_energy.o\
		quantum_energy.o\
		quantum_properties.o\
		velocity.o\
		hack_gms.o\
		hack_chrg.o\
		hack_vrc.o\
		hack_zbnd.o\
		energy.o\
		numerical_derivatives.o\
		dcd_io.o\
		multipoles.o\
		normal_mode_utilities.o\
		optimize_bfgs.o\
		optimize_coordinates.o\
		reaction_path.o\
		self_avoiding_walk.o\
		lbfgs.o\
		neb_climb.o\
		neb_spline.o\
		dcd_analysis.o\
		normal_mode.o\
		thermodynamics_RRHO.o\
		hack_urc.o\
		hack_gh.o\
		dynamics_utilities.o\
		dynamics_langevin_verlet.o\
		dynamics_leapfrog_verlet.o\
		dynamics_velocity_verlet.o\
		monte_carlo_energy.o\
		monte_carlo.o\
		quickdock.o\
		amberfile_io.o\
		my_random.o\
		lbfgsb.o cg+.o\
		lbfgsb_w.o\
		mini_wrapper.o\
		hack_kie.o\
		dynamo.o

.PHONY: all clean
.SUFFIXES: .F90 .F .c

#-------------------------------------------------------------------------------

all: $(OBJS)
	rm -f dynamo.a
	ar rUuc dynamo.a $(OBJS)

clean:
	rm -f *.o *.mod *.a

#-------------------------------------------------------------------------------

.F.o:
	$(FC) -c $(FLGS) $< -o $@

.F90.o:
	$(FC) -c $(FLGS) $< -o $@

.c.o:
	$(CC) -c $< -o $@
