# README

The function `plot_pca()` takes several arguments:
* __count_csv_file:__ A csv file with the raw_counts (not scaled) coming from an RNA-seq experiment. 
* __xp_design_csv_file:__ A file that links samples to their experimental condition.
* __timepoint:__ An integer for time (in hours) as indicated in the xp_design_csv_file. This will be used to filter the counts and keep only one timepoint.
* __pc_x_axis:__ = This has to be in the form "PCX" where X is the PC number. For instance, PC1 refers to the first component. This PC will be used to build the x axis. 
* __pc_y_axis__ = Similar to __pc_x_axis__ but this one will be used on the Y axis. 
* __pca_colour__ = An experimental factor present in the xp_design file. For instance, if "sequencing_batch" then the corresponding column will be used to label the points. 

The function 'do_permanova()' calculates the P value of your treatment effect. It takes several arguments:
* __count_csv_file:__ a csv file with the raw_counts (not scaled) coming from an RNA-seq experiment.
* __xp_design_csv_file:__ A file that links samples to their experimental condition (= treatment).
* __timepoint:__ An integer for time (in hours) as indicated in the xp_design_csv_file. This will be used to filter the counts and keep only one timepoint.
