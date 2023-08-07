"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import os

import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment

import hoomd_polymers
from hoomd_polymers.base.system import Pack
from hoomd_polymers.base.simulation import Simulation
from hoomd_polymers.library.polymers import LJChain
from hoomd_polymers.library.ff_classes import BeadSpring
import unyt


class MyProject(FlowProject):
    pass


class Borah(DefaultSlurmEnvironment):
    hostname_pattern = "borah"
    template = "borah.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="shortgpu",
            help="Specify the partition to submit to."
        )


class R2(DefaultSlurmEnvironment):
    hostname_pattern = "r2"
    template = "r2.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="shortgpuq",
            help="Specify the partition to submit to."
        )


class Fry(DefaultSlurmEnvironment):
    hostname_pattern = "fry"
    template = "fry.sh"

    @classmethod
    def add_args(cls, parser):
        parser.add_argument(
            "--partition",
            default="batch",
            help="Specify the partition to submit to."
        )

# Definition of project-related labels (classification)
@MyProject.label
def sim_done(job):
    return job.doc.sim_done


@MyProject.label
def sample_done(job):
    return job.doc.sample_done


def build_system(job):
    lj_chains = LJChain(
            num_mols=job.sp.num_mols,
            lengths=job.sp.chain_lengths,
            bead_sequence=job.sp.bead_sequence,
            bead_mass=job.sp.bead_mass,
            bond_lengths=job.sp.bond_lengths
    )

    system = Pack(
            molecules=lj_chains,
            density=job.sp.density,
            r_cut=job.sp.r_cut,
            packing_expand_factor=job.sp.packing_expand_factor,
            edge=job.sp.packing_edge,
    )

    length_units = getattr(unyt, job.sp.ref_length["units"])
    energy_units = getattr(unyt, job.sp.ref_energy["units"])
    mass_units = getattr(unyt, job.sp.ref_mass["units"])
    ref_length = job.sp.ref_length["value"] * length_units
    ref_energy = job.sp.ref_length["value"] * energy_units
    ref_mass = job.sp.ref_length["value"] * mass_units

    job.doc.ref_length = ref_length
    job.doc.ref_energy = ref_energy
    job.doc.ref_mass = ref_mass
    system.reference_length = ref_length
    system.reference_energy = ref_energy
    system.reference_mass = ref_mass
    print("------------")
    print(system._reference_values)
    print("------------")

	# Set up Forcefield:
    beads = {"A": dict(epsilon=job.sp.epsilon, sigma=job.sp.sigma)}
    bonds = {"A-A": dict(r0=job.sp.bond_r0, k=job.sp.bond_k)}
    angles = {"A-A-A": dict(t0=job.sp.angle_t0, k=job.sp.angle_k)}
    dihedrals = {
            "A-A-A-A": dict(
                phi0=job.sp.dihedral_phi0,
                k=job.sp.dihedral_k,
                d=-1,
                n=1
            )
    }

    bead_spring_ff = BeadSpring(
            beads=beads,
            bonds=bonds,
            angles=angles,
            dihedrals=dihedrals,
            r_cut=job.sp.r_cut
    )

    hoomd_ff = bead_spring_ff.hoomd_forcefield
    snapshot = system.hoomd_snapshot
    ref_units = (
            system.reference_length,
            system.reference_energy,
            system.reference_mass
    )
    #return snapshot, hoomd_ff, ref_units
    return system, hoomd_ff, ref_units


@MyProject.post(sim_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="simulation"
)
def run_sim(job):
    with job:
        # Create initial config, set up simulation
        system, hoomd_ff, ref_units = build_system(job)
        gsd_path = job.fn("trajectory.gsd")
        log_path = job.fn("log.txt")
        sim = Simulation(
                initial_state=system.hoomd_snapshot,
                forcefield=hoomd_ff,
                dt=job.sp.dt,
                r_cut=job.sp.r_cut,
                gsd_write_freq=job.sp.gsd_write_freq,
                gsd_file_name=gsd_path,
                log_write_freq=job.sp.log_write_freq,
                log_file_name=log_path,
        )
        sim.add_walls(wall_axis=(1,0,0), sigma=1, epsilon=1, r_cut=1.2)
        sim.pickle_forcefield(job.fn("forcefield.pickle"))
        sim.reference_length = ref_units[0]
        sim.reference_energy =  ref_units[1]
        sim.reference_mass = ref_units[2]
        tau_kT = sim.dt * job.sp.tau_kT
        target_box = (system.target_box /
                system.reference_length.value
        )
        # Store other sim information in job doc
        job.doc.wall_axis = (1,0,0)
        job.doc.wall_sigma = 1.0
        job.doc.wall_epsilon= 1.0
        job.doc.wall_r_cut = 1.2
        job.doc.tau_kT = tau_kT
        job.doc.target_box = target_box
        # Add time related job doc info
        job.doc.real_time_step = sim.real_timestep.to("fs").value
        job.doc.real_time_units = "fs"
        job.doc.n_steps = job.sp.n_steps
        job.doc.simulation_time = job.doc.real_time_step * job.doc.n_steps
        job.doc.gsd_step_time = job.doc.real_time_step * job.sp.gsd_write_freq
        job.doc.log_step_time = job.doc.real_time_step * job.sp.log_write_freq
        job.doc.shrink_time = job.doc.real_time_step * job.sp.shrink_n_steps
        # Set up stuff for shrinking volume step
        sim.save_restart_gsd(job.fn("init.gsd"))
        if job.sp.shrink_n_steps > 0:
            print("Running shrink step.")
            shrink_kT_ramp = sim.temperature_ramp(
                    n_steps=job.sp.shrink_n_steps,
                    kT_start=job.sp.shrink_kT,
                    kT_final=job.sp.kT
            )
            sim.run_update_volume(
                    final_box_lengths=target_box,
                    n_steps=job.sp.shrink_n_steps,
                    period=job.sp.shrink_period,
                    tau_kt=tau_kT,
                    kT=shrink_kT_ramp
            )
            print("Shrink step finished.")
        print("Running simulation.")
        sim.run_NVT(kT=job.sp.kT, n_steps=job.sp.n_steps, tau_kt=tau_kT)
        sim.save_restart_gsd(job.fn("restart.gsd"))
        print("Simulation finished.")
        job.doc.sim_done = True


@MyProject.pre(sim_done)
@MyProject.post(sample_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="sample"
)
def sample(job):
    from cmeutils.polymers import persistence_length
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("JOB TYPE: SAMPLE")
        print("------------------------------------")
        select_atoms_arg = " ".join(job.sp.bead_sequence)
        lp_mean, lp_std = persistence_length(
                job.fn("trajectory.gsd"),
                select_atoms_arg="name " + select_atoms_arg,
                window_size=50,
                start=-1000,
                stop=-1
        )
        job.doc.lp_mean = lp_mean
        job.doc.lp_std = lp_std
        print("Sampling finished.")
        job.doc.sampling_done = True


if __name__ == "__main__":
    MyProject().main()
