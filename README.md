# Pfizer_BioNTech_BNT162b2_vaccine_effectiveness_vs_Omicron
Code and data for "Predicted Symptomatic Effectiveness of  Pfizer-BioNTech BNT162b2 Vaccine  Against  Omicron Variant of SARS-CoV-2" paper 

Runs from the main.R file. 
Change "data_path <- 'YOUR_PATH/datafiles/'" in selectors.R as necessary. 

The code was developed and tested on MacOS in Pycharm R-plugin. 
The function showInChrome() in plotters.R shows plots as tabs in a Chrome browser window. It is  Mac and Chrome specific. 
If necessary, change this function or comment it out in plotters.R, but uncomment # fig, to show figures. 

Data sources: 
Crick_Comirnaty.csv 
https://github.com/davidlvb/Crick-UCLH-Legacy-VOCs-2021-05
Wall EC, Wu M, Harvey R, Kelly G, Warchal S, Sawyer C, Daniels R, Hobson P, Hatipoglu E, Ngai Y, Hussain S, Nicod J, Goldstone R, Ambrose K, Hindmarsh S, Beale R, Riddell A, Gamblin S, Howell M, Kassiotis G, Libri V, Williams B, Swanton C, Gandhi S, Bauer DLV. Neutralising antibody activity against SARS-CoV-2 VOCs B.1.617.2 and B.1.351 by BNT162b2 vaccination. Lancet. 2021 Jun 19;397(10292):2331-2333. doi: 10.1016/S0140-6736(21)01290-3. Epub 2021 Jun 3. PMID: 34090624; PMCID: PMC8175044.


Polymutant.csv 
https://www.nature.com/articles/s41586-021-04005-0
Schmidt, F., Weisblum, Y., Rutkowska, M. et al. High genetic barrier to SARS-CoV-2 polyclonal neutralizing antibody escape. Nature (2021). https://doi.org/10.1038/s41586-021-04005-0
