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
def seed_done(job):
    return job.isfile("seed-restart.gsd")


@MyProject.label
def sample_done(job):
    return job.doc.sample_done


@MyProject.label
def seed_system(job):
    return job.sp.n_duplicates is None


@MyProject.label
def not_seed_system(job):
    return job.sp.n_duplicates is not None

def duplicate_seed(job):
    from flowermd.modules.welding import Interface
    n_dups = 0
    while n_dups != job.sp.n_duplicates:
        for i, axis in enumerate([(1,0,0), (0,1,0), (0,0,1)]):
            if n_dups == job.sp.n_duplicates:
                break
            dup = Interface(
                    gsd_files=[job.fn("seed-restart.gsd")],
                    interface_axis=axis,
                    gap=8,
                    remove_void_particles=False
            )
            with gsd.hoomd.open(job.fn("seed-restart.gsd"), "w") as traj:
                traj.append(dup.hoomd_snapshot)
            n_dups += 1

    with gsd.hoomd.open(job.fn("seed-restart.gsd"), "r") as traj:
        job.doc.N = traj[0].particles.N


@MyProject.pre(seed_system)
@MyProject.post(seed_done)
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="nvt"
)
def make_seed(job):
    import unyt as u
    from unyt import Unit
    import flowermd
    from flowermd.base.system import Pack
    from flowermd.base.simulation import Simulation
    from flowermd.library import PEKK_Para, GAFF 
    from flowermd.utils import get_target_box_mass_density

    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
    
        molecules = PEKK_Para(lengths=job.sp.lengths, num_mols=job.sp.num_mols)
        system = Pack(molecules=molecules, density=job.sp.density) 

        system.apply_forcefield(
                force_field=OPLS_AA(),
                r_cut=job.sp.r_cut,
                auto_scale=job.sp.auto_scale,
                scale_charges=True,
                remove_charges=job.sp.remove_charges,
                remove_hydrogens=job.sp.remove_hydrogens,
                pppm_resolution=job.sp.pppm_resolution,
                pppm_order=job.sp.pppm_order
        )

        gsd_path = job.fn("seed.gsd")
        log_path = job.fn("seed-log.txt")

        sim = Simulation.from_system(
                system=system,
                dt=job.sp.dt,
                gsd_write_freq=job.sp.gsd_write_freq,
                gsd_file_name=gsd_path,
                log_write_freq=job.sp.log_write_freq,
                log_file_name=log_path,
        )
        sim.pickle_forcefield(job.fn("forcefield.pickle"))
        # Store unit information in job doc
        tau_kT = sim.dt * job.sp.tau_kT
        job.doc.tau_kT = tau_kT
        job.doc.ref_mass = sim.reference_mass.to("amu").value
        job.doc.ref_mass_units = "amu"
        job.doc.ref_energy = sim.reference_energy.to("kJ/mol").value
        job.doc.ref_energy_units = "kJ/mol"
        job.doc.ref_length = sim.reference_length.to("nm").value
        job.doc.ref_length_units = "nm"
        job.doc.real_time_step = sim.real_timestep.to("fs").value
        job.doc.real_time_units = "fs"
        # Set up stuff for shrinking volume step
        print("Running shrink step.")
        shrink_kT_ramp = sim.temperature_ramp(
                n_steps=job.sp.shrink_n_steps,
                kT_start=job.sp.shrink_kT,
                kT_final=job.sp.kT
        )
        target_box = get_target_box_mass_density(
                mass=system.mass.to("g"),
                density=job.sp.density * (Unit("g/cm**3") / 2)
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
        sim.run_NVT(kT=6.0, n_steps=2e6, tau_kt=tau_kT)
        sim.save_restart_gsd(job.fn("seed-restart.gsd"))
        print("Seed simulation finished.")
        duplicate_seed(job)
        print("Seed duplication finished.")


@MyProject.pre(nvt_done)
@MyProject.post(sample_done)
@MyProject.operation(
        directives={"ngpu": 0, "executable": "python -u"}, name="sample"
)
def sample(job):
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here
        job.doc.sample_done = True


if __name__ == "__main__":
    MyProject(environment=Fry).main()
