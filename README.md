The following MATLAB 2024 scripts are provided:

**incremental_stress** - A self-contained stress-prescribing incremental simulation script. 
                         This can be used to simulate multiaxial and non-proportional loadings.

**incremental_strain** - A self-contained strain-prescribing incremental simulation script. 
                         This script effectively guesses the stresses that will lead to the specified plastic strains. 
                         Simplifying assumptions of zero hydrostatic stress and stress-increments normal to the yield 
                         surface have been made. Again, multiaxial and non-proportional loadings can be simulated.

**PLOT_pi_plane**      - A script constructed to visualise the results of the above two simulation scripts.
                         The plastic development is plotted in the pi-plane, and GIF's can be constructed.

**total_form**         - The script used in parameter optimisations. It is a function, and cannot be called by itself.

**runparams**          - An initialisation and visualisation script made to run the total_form function based on
                         experimental data.

**NonlinLeastSquares** - A self-contained script used to optimise isotropic models like the Voce and Hollomon.

**MultiFullFit**       - The optimisation script used to run fmincon() for Chaboche parameter fits. 
                         The script calls the total_form function around 1000 times during an optimisation.

Some example material data has been provided, as running the scripts is otherwise difficult.
The data has been modified, and does therefore not represent real material behaviour.
