V.T. Cooper, Univ of Washington, 29 Sep 2021
-----------
Analysis supporting Cooper, Roach, Thomson, Brenner, Smith, and Bitz (2021).

Processes hourly coupling version of Roach et al. 2019 model.
NSIDC/NOAA Sea Ice Concentration and in situ observations from BGOS and SODA, 
supplemented by Arctic Sea State SWIFT buoys, all from the Beaufort Sea.

For the paper, the jgr_*.ipynb notebooks have the preprocessing and analysis 
to support all figures in paper.

Chronology of scratch analysis:
waveice_analysis.ipynb, 
then wave_ice_analysis_2021_pub.ipynb are the main notebooks. The latter
was cleaned up to become jgr_waveice_analysis.ipynb.

waveice_climatology.ipynb was used to output climatology of model results
and waveice_coupled_2021.ipynb was used to preprocess the new hourly run.

model_dev_analysis.ipynb is separately used to explore new experiments.

The convert-mat* files were used to preprocess the original .mat data from
APL team.

Some misc. interim outputs are saved here to save time when running analysis 
notebook.

