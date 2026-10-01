config = {
    "problem": "CONTACT_move",
    "eos_list": [
        {
            "eos": "ideal",
            "gamma": 1.4,
        },
        {
            "eos": "ideal",
            "gamma": 1.4,
        },
    ],
    "kappa": 1.0e-100,
    "mean_molecular_weight": 1.0,
    "problem_parameters": {
        "convergence_order": 1,
        "convergence_metric": "full",
        "exclusion_width": 0.1,
    },
    "t0": 0.0,
    "tf": 0.1,
    "max_steps": 15000,
    "cfl": 0.95,
    "max_dt": 0.01,
    "log_frequency": 10,
    "output_frequency": 1000,
    "levels": [9, 1, 1],
    "coords": [
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
    ],
    "boundaries": {
        "xlow": "flow",
        "xhigh": "flow",
        "ylow": "flow",
        "yhigh": "flow",
        "zlow": "flow",
        "zhigh": "flow",
    },
    "output_method": "csv",
}