#!/usr/bin/env python
"""Initialize the project's data space.

Iterates over all defined state points and initializes
the associated job workspace directories.
The result of running this file is the creation of a signac workspace:
    - signac.rc file containing the project name
    - signac_statepoints.json summary for the entire workspace
    - workspace/ directory that contains a sub-directory of every individual statepoint
    - signac_statepoints.json within each individual statepoint sub-directory.

"""

import signac
import flow
import logging
from collections import OrderedDict
from itertools import product


def get_parameters():
    ''''''
    # System Details:
    parameters = OrderedDict()
    parameters["bead_sequence"] = [["A"]]
    parameters["chain_lengths"] = [60]
    parameters["num_mols"] = [1]
    parameters["bead_mass"] = [{"A": 450}]
    parameters["bond_lengths"] = [{"A-A": 1.0}]
    parameters["density"] = [0.01]
    parameters["ref_length"] = [dict(value=1, units="nm")]
    parameters["ref_energy"] = [dict(value=1, units="kJ")]
    parameters["ref_mass"] = [dict(value=450, units="amu")]
    parameters["packing_expand_factor"] = [8]
    parameters["packing_edge"] = [0.5]
    # Forcefield Details:
    parameters["sigma"] = [1.0]
    parameters["epsilon"] = [1.0]
    # Bonds
    parameters["bond_k"] = [500]
    parameters["bond_r0"] = [1.1]
    # Angles
    parameters["angle_k"] = [200, 300]
    parameters["angle_t0"] = [1.6, 1.8, 2.0, 2.2, 2.4, 2.6, 2.8, 3.0]
    # Dihedrals
    parameters["dihedral_k"] = [5, 50, 100]
    parameters["dihedral_phi0"] = [0.0]
    # Other
    parameters["r_cut"] = [2.5]
    parameters["exclusions"] = [["bond", "1-3"]]

    # Simulation Details:
    parameters["kT"] = [8.0]
    parameters["n_steps"] = [2e7]
    parameters["shrink_kT"] = [8.0]
    parameters["shrink_n_steps"] = [0]
    parameters["shrink_period"] = [1]
    parameters["dt"] = [0.0003]
    parameters["tau_kT"] = [100]
    parameters["gsd_write_freq"] = [1e4]
    parameters["log_write_freq"] = [1e3]
    return list(parameters.keys()), list(product(*parameters.values()))


def main():
    project = signac.init_project() # Set the signac project name
    param_names, param_combinations = get_parameters()
    # Create the generate jobs
    for params in param_combinations:
        statepoint = dict(zip(param_names, params))
        job = project.open_job(statepoint)
        job.init()
        job.doc.setdefault("sim_done", False)
        job.doc.setdefault("sample_done", False)
    print(f"Creating {len(param_combinations)} job workspaces...")


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    main()
