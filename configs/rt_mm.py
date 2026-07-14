config = {
    "problem": "rt-test_mm",
    "eos_list": [
        { # bottom lighter fluid
            "eos": "ideal",
            "gamma": 1.4,
        },
        { # top heavier fluid
            "eos": "ideal",
            "gamma": 1.4,
        },
    ],
    "kappa": 0.0,
    "adaptive_check": False,
    "limiter_id": 0,
    "closure_id": 0,
    "gravity_acc": [0.0, -0.1, 0.0],
    "mean_molecular_weight": 1.0,
    "t0": 0.0,
    "tf": 50.0,
    "max_steps": 100000,
    "cfl": 0.75,
    "max_dt": 5.0e-2,
    "log_frequency": 1000,
    "output_frequency": 1000,
    "levels": [7, 9, 4],
    "coords": [
        [0.0, 0.0, 0.0],  # [-0.25, -0.75, 0.0],
        [0.5, 1.4, 1.0],  # [0.25, 0.75, 1.0],
    ],
    "boundaries": {
        "xlow": "periodic",  # needs to be changed to Homogeneous Neumann
        "xhigh": "periodic",  # needs to be changed to Homogeneous Neumann
        "ylow": "reflecting",
        "yhigh": "reflecting",
        "zlow": "flow",
        "zhigh": "flow",
    },
    "output_method": "vti",
}
