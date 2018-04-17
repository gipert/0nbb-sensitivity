# 0nbb-sensitivity
The `discSensitivityVsBI` executable calculates the discovery sensitivity (half-life for which one correctly reconstructs a signal in 50% of simulated experiments) to a gaussian signal over a flat background vs. the background index `BI`.

## Usage tips
* Compile everything with `make -j4`. Thread parallelization with OpenMP is currently enabled only on Linux systems.
* Control the number of threads with the environment variables `OMP_NUM_THREADS` and `CUBACORES`, respectively for the number of threads of the main loop and for the number of threads used by Cuba during the integration.
* All the settings you may need to adjust are stored in `masterconf.json`
