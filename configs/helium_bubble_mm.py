config = {
    "problem": "helium_bubble_mm",

    # Keep ideal EOS for both materials.
    # material 0 = ambient / shock-tube gas
    # material 1 = helium bubble
    "eos_list": [
        {
            "eos": "ideal",
            "gamma": 1.4,
        },
        {
            "eos": "ideal",
            "gamma": 1.648,
        },
    ],

    "kappa": 0.0,
    "mean_molecular_weight": 1.0,

    "t0": 0.0,
    "tf": 2.5,
    "max_steps": 20000,

    "cfl": 0.5,
    "max_dt": 0.01,

    "log_frequency": 20,
    "output_frequency": 20,

    # Previous large physical domain:
    # "levels": [8, 6, 4],
    # "coords": [
    #     [0.0, -40.0, 0.0],
    #     [350.0, 40.0, 1.0],
    # ],
    #
    # New simple unit-box problem.
    "levels": [9, 8, 6],
    "lowest_level": 6,

    "coords": [
        [0.0, 0.0, 0.0],
        [4.0, 2.0, 2.0],
    ],

    "boundaries": {
        "xlow": "flow",
        "xhigh": "flow",
        "ylow": "reflecting",
        "yhigh": "reflecting",
        "zlow": "reflecting",
        "zhigh": "reflecting",
    },

    "output_method": "vti",
}
