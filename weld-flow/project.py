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

def get_slab_file(job, file_name):

    os.path.join(job.path, "..", "..", "slab-files"


@MyProject.post(sim_done)
<<<<<<<< HEAD:weld-flow/project.py
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="weld"
)
def run_weld(job):
|||||||| 598f896:src/project.py
@MyProject.operation(directives={"ngpu": 1, "executable": "python -u"})
def run_simulation(job):
========
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="simulate"
)
def run_simulation(job):
>>>>>>>> main:flow-templates/basic/project.py
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here


@MyProject.pre(sim_done)
@MyProject.post(sample_done)
<<<<<<<< HEAD:weld-flow/project.py
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="weld"
)
|||||||| 598f896:src/project.py
@MyProject.operation(directives={"ngpu": 1, "executable": "python -u"})
========
@MyProject.operation(
        directives={"ngpu": 1, "executable": "python -u"}, name="sample"
)
>>>>>>>> main:flow-templates/basic/project.py
def sample(job):
    # Add package imports here
    with job:
        print("JOB ID NUMBER:")
        print(job.id)
        print("------------------------------------")
        # Add your script here


if __name__ == "__main__":
    MyProject().main()
