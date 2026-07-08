config = {
    "problem": "kh-test_mm",
    "eos_list": [
        {
            "eos": "ideal",
            "gamma": 1.6667,
        },
        {
            "eos": "ideal",
            "gamma": 1.6667,
        },
        ],
    "kappa": 0.0,
    "mean_molecular_weight": 1.0,
    "eos": "ideal",
    "t0": 0.0,
    "tf": 1.0,
    "max_steps": 10000,
    "cfl": 0.5,
    "max_dt": 0.001,
    "log_frequency": 100,
    "output_frequency": 100,
    "levels": [8, 8, 4],
    "lowest_level": 6,
    "coords": [
        [0.0, 0.0, 0.0],
        [1.0, 1.0, 1.0],
    ],
    "boundaries": {
        "xlow": "periodic",
        "xhigh": "periodic",
        "ylow": "periodic",
        "yhigh": "periodic",
        "zlow": "flow",
        "zhigh": "flow",
    },
    "output_method": "vti",
}
