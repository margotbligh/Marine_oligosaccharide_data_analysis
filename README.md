# Marine_oligosaccharide_data_analysis
Code used to analyse LC-HRMS data described in marine oligosaccharide paper

The paper: insert link and citation

LC-HRMS data (mzML) format is available on the MassIVE data base. 

Additional data objects:

- my_own_model_10000-lowlr_016.h5 = model trained for the study
- R data objects from pre-processing (sequential)
  - data.peaks.RData = R data object with picked peaks (output of findChromPeaks)
  - data.peaks.ref.RData = R data object with refined peaks (output of refineChromPeaks)
  - data.peaks.ref.gsub.RData = R data object with initial peak grouping for subset based retention time alignment
  - data.peaks.ref.gsub.rt.RData = R data object with aligned peaks (output of adjustRtime)
  - data.peaks.ref.gsub.rt.grp.RData = R data object with grouped peaks (output of groupChromPeaks)
  - data.peaks.ref.gsub.rt.grp.fill.RData = R data object with filled peaks (output of fillChromPeaks)
- tables
  - aligned_feature_table_deiso.csv = feature table in NeatMS input format after deisotoping - **used as input for NeatMS denoising**
  - neatms_export_with_extra_properties.csv = OUTPUT of NeatMS denoising
  - peaks_manually_checked.csv = list of retained peaks after manual review stage



