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
            "gamma": 1.67,
        },
    ],

    "kappa": 0.0,
    "mean_molecular_weight": 1.0,

    "t0": 0.0,
    "tf": 6,
    "max_steps": 30000,

    "cfl": 0.8,
    "max_dt": 0.01,

    "log_frequency": 50,
    "output_frequency": 100,

    "levels": [9, 8, 8],
    "lowest_level": 8,

    "coords": [
        [0.0, 0.0, 0.0],
        [352.0, 89.0, 89.0],
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
    "extp_steps": 30, # used for extrapolation
    "debug_extrapolation_t0": False,
    "gfm_method": "riemann_gfm",
}
