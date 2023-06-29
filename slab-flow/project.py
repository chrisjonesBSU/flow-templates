"""Define the project's workflow logic and operation functions.

Execute this script directly from the command line, to view your project's
status, execute operations and submit them to a cluster. See also:

    $ python src/project.py --help
"""
import signac
from flow import FlowProject, directives
from flow.environment import DefaultSlurmEnvironment
import os


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


@MyProject.post(sim_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="sim"
)
def run_sim(job):
    import hoomd_polymers
    from hoomd_polymers.base.system import Pack
    from hoomd_polymers.base.simulation import Simulaton 
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        mol_obj = getattr(hoomd_polymers.library.polymers, job.sp.molecule)
        mols = mol_obj(
                num_moles=job.sp.num_mols,
                lengths=job.sp.lengths,
                force_field=job.sp.forcefield
        )
        system = Pack(
                molecules=[mols],
                density=job.sp.density,
                r_cut=2.5,
                auto_scale=True
        ) 

        if job.sp.remove_hydrogens:
            system.remove_hydrogens()
        if job.sp.remove_charges:
            system.remove_charges()

        sim = Simulation(
                initial_state=system.hoomd_snapshot,
                forcefield=system.hoomd_forcefield,
                dt=job.sp.dt,
                gsd_write_freq=job.sp.gsd_write_freq,
                log_write_freq=job.sp.log_write_freq,
        )
        sim.reference_length = system.reference_length
        sim.reference_energy = system.reference_energy
        sim.reference_mass = system.reference_mass
        # Store unit information in job doc
        tau_kT = sim.dt * job.sp.tau_kT 
        job.doc.tau_kT = tau_kT
        job.doc.target_box = target_box
        job.doc.ref_mass = sim.reference_mass.to("amu").value()
        job.doc.ref_mass_units = "amu"
        job.doc.ref_energy = sim.reference_energy.to("kJ/mol").value()
        job.doc.ref_energy_units = "kJ/mol"
        job.doc.ref_length = sim.reference_length.to("angstrom").value()
        job.doc.ref_length_units = "angstrom"
        job.doc.real_time_step = sim.real_timestep.to("fs").value()
        job.doc.real_time_units = "fs"
        job.doc.simulation_time = job.doc.real_time_step * job.sp.n_steps
        job.doc.shrink_time = job.doc.real_time_step * job.sp.shrink_n_steps
        # Add wall forces
        sim.add_walls(wall_axis=(1,0,0), sigma=1.0, epsilon=1.0, r_cut=2.5)
        # Set up stuff for shrinking volume step 
        print("Running shrink step...")
        target_box = system.target_box * system.reference_distance.to("angstrom").value()
        shrink_kT_ramp = sim.kT_ramp(
                n_steps=job.sp.shrink_n_steps,
                kT_start=job.sp.shrink_kT,
                kT_final=job.sp.kT
        )
        sim.run_update_volume(
                final_box=target_box,
                n_steps=job.sp.shrink_n_steps,
                period=job.sp.shrink_period,
                tau_kt=tau_kT,
                kT=shrink_kT_ramp
        )
        print("Shrink step finished...")
        print("Running simulation")
        sim.run_NVT(kT=job.sp.kT, n_steps=job.sp.n_steps, tau_kt=tau_kT)
                



@MyProject.pre(sim_done)
@MyProject.post(sample_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="weld"
)
def sample(job):
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here


if __name__ == "__main__":
    MyProject().main()
