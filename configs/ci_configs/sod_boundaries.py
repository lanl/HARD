config = {
      "problem": "sod",
      "gamma": 1.4,
      "eos": "ideal",
      "kappa": 1.0e-100,
      "mean_molecular_weight": 1.0,
      "t0": 0.0,
      "tf": 0.07,
      "max_steps": 500,
      "cfl": 0.3,
      "max_dt": 1.0,
      "log_frequency": 500,
      "output_frequency": 500,
      "levels": [10, 3, 3],
      "coords": [
          [0.0, 0.0, 0.0],
          [1.0, 1.0, 1.0],
      ],
      # Per-field boundaries: the left and right Sod states are pinned with
      # dirichlet values, the vectors take reflecting on the right. Identical
      # to plain flow at tf=0.07 since no wave reaches a boundary.
      "boundaries": {
          "xlow": [
              "flow",
              ("mass_density", "dirichlet", 1.0),
              ("momentum_density", "dirichlet", [0.0, 0.0, 0.0]),
              ("total_energy_density", "dirichlet", 2.5),
              ("velocity", "dirichlet", [0.0, 0.0, 0.0]),
              ("pressure", "dirichlet", 1.0),
              ("specific_internal_energy", "dirichlet", 2.5),
          ],
          "xhigh": [
              "flow",
              ("mass_density", "dirichlet", 0.125),
              ("momentum_density", "reflecting"),
              ("total_energy_density", "dirichlet", 0.25),
              ("velocity", "reflecting"),
              ("pressure", "dirichlet", 0.1),
              ("specific_internal_energy", "dirichlet", 2.0),
          ],
          "ylow": "flow",
          "yhigh": "flow",
          "zlow": "flow",
          "zhigh": "flow",
      },
      "output_method": "csv",
  }
