config = {
    "problem": "sod_mm",

    "eos_list": [
        {
            "eos": "ideal",
            "gamma": 1.4,
        },
        {
            "eos": "ideal",
            "gamma": 2.0,
        },
    ],
    # # --------------------------------------------------- # here 2 materials
    # # add here the things for the different materials
    # # --------------------------------------------------- #
    "kappa": 1.0e-100,
    "mean_molecular_weight": 1.0,
    "t0": 0.0,
    "tf": 0.3,
    "max_steps": 200,
    "cfl": 0.3,
    "max_dt": 0.0001,
    "log_frequency": 50,
    "output_frequency": 500,
    "levels": [8, 8, 8],
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
    "extp_steps": 30, # used for extrapolation
    "output_method": "csv",
}
