from collections import namedtuple
from typing import *
import numpy as np
import mdsine2 as md2


GLVParamSet = namedtuple(
    'GLVParamSet',
    [
        'growth', 'interactions', 'perturbations'
    ]
)


def generate_mouse_name(mouse_index: int) -> str:
    return "SYNTH_SUBJ_{}".format(mouse_index)


def subj_with_most_timepoints(study: md2.Study) -> md2.Subject:
    best_subj = None
    for subj in study:
        if best_subj is None:
            best_subj = subj
        elif len(best_subj.times) < len(subj.times):
            best_subj = subj
    if best_subj is None:
        raise ValueError("No subjects found in study.")
    return best_subj


def forward_simulate_subject(
        glv_params: GLVParamSet,
        pert_starts: List[float],
        pert_ends: List[float],
        time_points: List[float],
        initial_conditions: np.ndarray,
        dt: float,
        sim_max: float,
        time_subset: bool = True
) -> Tuple[np.ndarray, np.ndarray]:
    """Forward simulation for a single subject"""
    dyn = md2.model.gLVDynamicsSingleClustering(
        growth=glv_params.growth,
        interactions=glv_params.interactions,
        perturbations=glv_params.perturbations,
        perturbation_starts=pert_starts,
        perturbation_ends=pert_ends,
        start_day=time_points[0],
        sim_max=sim_max
    )

    if len(initial_conditions.shape) == 1:
        initial_conditions = np.expand_dims(initial_conditions, 1)
    if time_subset:
        x = md2.integrate(
            dynamics=dyn,
            initial_conditions=initial_conditions,
            dt=dt,
            final_day=time_points[-1],
            subsample=True,
            times=time_points
        )
    else:
        x = md2.integrate(
            dynamics=dyn,
            initial_conditions=initial_conditions,
            dt=dt,
            final_day=time_points[-1],
            subsample=False
        )
    fwsim_values = x['X']
    times = x['times']
    return fwsim_values, times
