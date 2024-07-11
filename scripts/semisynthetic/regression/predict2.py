import argparse
from pathlib import Path

import numpy as np
import pickle
from scipy.integrate import solve_ivp
import mdsine2 as md2
from numba import jit
from tqdm import tqdm


def forward_sim_single_subj_glv(A, g, B, x0, u, times, rel_abund=False):
    """
    forward simulate for a single subject
    (np.ndarray) x0 : N dimensional array containing the initial abundances(log)
    (np.ndarray) A, g, B: : the coefficients for interactions, growth, perturbation
    (np.ndarray) u : the perturbation indicator
    (np.ndarray) t : the time coefficients
    """
    # def grad_fn(A, g, B, u):
    #     def fn(t, x):
    #         if B is None or u is None:
    #             return g + A.dot(x)
    #         elif B is not None and u is not None:
    #             return g + A.dot(np.exp(x)) + B.dot(u)
    #     return fn

    def grad_fn(A, g, B, u):
        def fn1(t, x):
            return g + A.dot(x)
        def fn2(t, x):
            return g + A.dot(np.exp(x)) + B.dot(u)
        if B is None or u is None:
            return jit(fn1)
        elif B is not None and u is not None:
            return jit(fn2)
        else:
            raise Exception("Unknown function class")

    x_pred = np.zeros((times.shape[0], x0.shape[0]))
    x_pred[0] = np.exp(x0)
    xt = x0
    if np.sum(B) == 0:
        B = None

    for t in tqdm(range(1, times.shape[0])):
        if u is not None:
            grad = grad_fn(A, g, B, u[t-1])
        else:
            grad = grad_fn(A, g, None, None)
        dt = times[t] - times[t-1]
        ivp = solve_ivp(grad, (0, 0+dt), xt, method="RK45")
        xt = ivp.y[:, -1]
        if rel_abund:
            x_pred[t] = np.exp(xt) / np.sum(np.exp(xt))
        else:
            x_pred[t] = np.exp(xt)

    return x_pred


def forward_sim_glv(data_path: Path,
                    x0: np.ndarray,
                    u: np.ndarray,
                    sim_times: np.ndarray) -> np.ndarray:
    print(f"Evaluating gLV simulation using output ({data_path})")
    with open(data_path, "rb") as f:
        model = pickle.load(f)
        A, g, B = model.get_params()

    x0 = np.copy(x0)
    """
    gLV inference is run with scaling (for numerical precision).
    This is the inverse transformation!
    """
    x0 = np.log(x0)

    # Include the limit of detection value
    return forward_sim_single_subj_glv(A, g, B, x0, u, sim_times, rel_abund=False).transpose(1, 0)


def perturbations_to_u(study: md2.Study, subject: md2.Subject, times: np.ndarray) -> np.ndarray:
    n_perts = len(study.perturbations)
    perturb_ids = np.zeros((len(times), n_perts), dtype=int)
    for p_idx, pert in enumerate(study.perturbations):
        start = pert.starts[subject.name]
        end = pert.ends[subject.name]
        timepoint_indices, = np.where((times >= start) & (times <= end))
        perturb_ids[timepoint_indices, p_idx] = 1
    return perturb_ids


def main(glv_pkl_path: Path, study: md2.Study, start_time: float, x0: np.ndarray, x0_scaling: float, final_day: float, sim_dt: float, out_path: Path):
    sim_times = np.arange(start=start_time, stop=final_day + sim_dt, step=sim_dt)
    u = perturbations_to_u(study, subject, sim_times)
    fwsims = forward_sim_glv(glv_pkl_path, x0 * x0_scaling, u, sim_times)
    np.savez(out_path, sims=fwsims / x0_scaling, times=sim_times)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument('--glv-pkl', dest='glv_pkl_path', type=str)
    parser.add_argument('--regression-inputs-dir', dest='regression_inputs_dir', type=str)
    parser.add_argument('--study', dest='study', type=str)
    parser.add_argument('--subject', dest='subject', type=str)
    parser.add_argument('--out-path', dest='out_path', type=str)
    parser.add_argument('--sim-dt', dest='sim_dt', type=float, default=0.01)
    parser.add_argument('--sim-max', dest='sim_max', type=float, default=1e20)
    parser.add_argument('--initialize_from_study', type=str, dest='init_study',
                        required=False,
                        help='[Experimental, for semisynthetic comparison only --  NOT TO BE COMMMITED TO GIT] The path to another study file which contains data different from the '
                             'specified study. This study\'s taxa must be a superset of the target study taxa.'
                        )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()

    study = md2.Study.load(args.study)
    subj_name = args.subject
    subject = None
    subject_idx = -1

    for subj_idx, subj in enumerate(study):
        if subj.name == subj_name:
            subject = subj
            subject_idx = subj_idx
            break
    if subject is None:
        print("Unknown subject `{}`. Available subjects: [{}]".format(
            subj_name,
            ",".join([subj.name for subj in study])
        ))
        exit(1)

    # ======= Initial conditions
    if args.init_study is None:
        M = subject.matrix()['abs']
        initial_conditions = M[:, 0]
    else:
        other_study = md2.Study.load(args.init_study)
        other_taxa = other_study.taxa
        M = other_study[subject.name].matrix()['abs']
        initial_conditions = M[:, 0]

        taxa_indices = [other_taxa[taxa.name].idx for taxa in study.taxa]
        initial_conditions = initial_conditions[taxa_indices] + 1e4

    abund_scaling = np.load(Path(args.regression_inputs_dir) / 'rescale.np')
    main(
        glv_pkl_path=Path(args.glv_pkl_path),
        study=study,
        out_path=Path(args.out_path),
        start_time=subject.times[0],
        x0=initial_conditions,
        x0_scaling=abund_scaling,
        final_day=subject.times[-1],
        sim_dt=args.sim_dt
    )