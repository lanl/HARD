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
    "t0": 0.0,
    "tf": 0.2,
    "max_steps": 10000,
    "cfl": 0.3,
    "max_dt": 0.01,
    "log_frequency": 500,
    "output_frequency": 1000,
    "levels": [9, 9, 9],
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
    "catalyst": {
        "script": "/path",
        "implementation": "paraview",
        "implementation_directory": "/path",
    },
}
