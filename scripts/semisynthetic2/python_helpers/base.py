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


# def forward_simulate(
#         glv_params: GLVParamSet,
#         study: md2.Study,
#         dt: float,
#         sim_max: float
# ) -> Dict[str, np.ndarray]:
#     """Forward simulation for all subjects in a Study object"""
#     return {
#         subj.name: forward_simulate_subject(glv_params, study, subj, dt, sim_max)[0]
#         for subj in study
#     }


def forward_simulate_subject(
        glv_params: GLVParamSet,
        pert_starts: List[float],
        pert_ends: List[float],
        time_points: List[float],
        initial_conditions: np.ndarray,
        dt: float,
        sim_max: float,
        time_subset: bool = True
) -> np.ndarray:
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

    if len(initial_conditions) == 1:
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
