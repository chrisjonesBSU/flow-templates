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
        directives={"ngpu": 1, "executable": "python -u"}, name="npt"
)
def run_npt(job):
    import hoomd_polymers
    from hoomd_polymers.base.system import Pack
    from hoomd_polymers.base.simulation import Simulaton 
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        mol_cls = getattr(hoomd_polymers.library.polymers, job.sp.molecule)
        ff = getattr(hoomd_polymers.forcefields, job.sp.forcefield)
        mol_obj = mol_cls(
                    num_moles=job.sp.num_mols,
                    lengths=job.sp.lengths,
                    force_field=ff()
                )

        system = Pack(
                    molecule=mol_obj,
                    density=job.sp.density,
                    r_cut=job.sp.r_cut,
                    auto_scale=True,
                    remove_hydrogens=job.sp.remove_hydrogens,
                    remove_charges=job.sp.remove_charges
                ) 

        gsd_path = job.fn("trajectory.gsd")
        log_path = job.fn("log.txt")
        sim = Simulation(
                initial_state=system.hoomd_snapshot,
                forcefield=system.hoomd_forcefield,
                dt=job.sp.dt,
                r_cut=job.sp.r_cut,
                gsd_write_freq=job.sp.gsd_write_freq,
                gsd_file_name=gsd_path,
                log_write_freq=job.sp.log_write_freq,
                log_file_name=log_path,
        )
        sim.pickle_forcefield(job.fn("forcefield.pickle"))
        sim.save_restart_gsd(job.fn("init.gsd"))
        sim.reference_length = system.reference_length * job.sp.sigma_scale
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
        job.doc.ref_length = sim.reference_length.to("nm").value()
        job.doc.ref_length_units = "nm"
        job.doc.real_time_step = sim.real_timestep.to("fs").value()
        job.doc.real_time_units = "fs"
        job.doc.n_steps = job.sp.n_steps
        job.doc.simulation_time = job.doc.real_time_step * job.doc.n_steps
        job.doc.shrink_time = job.doc.real_time_step * job.sp.shrink_n_steps
        # Set up stuff for shrinking volume step 
        print("Running shrink step.")
        target_box = system.target_box/system.reference_distance.value()
        shrink_kT_ramp = sim.kT_ramp(
                n_steps=job.sp.shrink_n_steps,
                kT_start=job.sp.shrink_kT,
                kT_final=job.sp.kT
        )

        # Anneal to just below target density
        sim.run_update_volume(
                final_box=target_box*0.85,
                n_steps=job.sp.shrink_n_steps,
                period=job.sp.shrink_period,
                tau_kt=tau_kT,
                kT=shrink_kT_ramp
        )

        # Run for a bit at lower density
        sim.run_NVT(n_steps=2e7, kT=job.sp.kT, tau_kt=tau_kT)

        # Compress
        sim.run_update_volume(
                final_box=target_box*1.15,
                n_steps=5e6,
                period=1000,
                tau_kt=tau_kT,
                kT=job.sp.kT
        )

        # Run for a bit at higher density
        sim.run_NVT(n_steps=2e7, kT=job.sp.kT, tau_kt=tau_kT)

        # Expand back to target density
        sim.run_update_volume(
                final_box=target_box,
                n_steps=5e6,
                period=1000,
                tau_kt=tau_kT,
                kT=job.sp.kT
        )
        # Short run at NVT
        sim.run_NVT(n_steps=2e7, kT=job.sp.kT, tau_kt=tau_kT)
        print("Shrink step finished.")
        print("Running NPT simulation.")
        sim.run_NPT(
                kT=job.sp.kT,
                pressure=job.sp.pressure,
                n_steps=job.sp.n_steps,
                tau_kt=tau_kt,
                tau_pressure=job.sp.tau_pressure*sim.dt
        )
        sim.save_restart_gsd(job.fn("npt-restart.gsd"))
        print("Simulation finished.")


@MyProject.pre(sim_done)
@MyProject.post(sample_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="sample-npt"
)
def sample_npt(job):
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here


if __name__ == "__main__":
    MyProject().main()
