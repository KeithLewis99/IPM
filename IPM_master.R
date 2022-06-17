
# The idea here is to have a loop for each model variation and then send the results to a new dashboard. I can get it to loop but not make different versions of the dashboard.  For right now, just change "b" here and in IPM_out and it creates the dashboard with all the figures in new folders.  

# IPM_out_diag----
b <- c(1:3)
b = 2

# folder files should follow this format: model #, time, description
for (i in b){
     outfile <- paste0("IPM_out_diag", i, ".html")
     outDir1 <- paste0("output/M29_1985P_20K_CPreg_b=", i)
     #outDir2 <- paste0("output/output_ez", i)
     #if(!dir.exists(paste0("output/", outDir)))dir.create(outDir)
     rmarkdown::render("IPM_out_diag.Rmd", output_file = outfile, output_dir = outDir1) # makes a dashboard
     #ezknitr::ezknit(file = "IPM_out_diag.Rmd", out_dir = outDir2, keep_md = F) # saves all graphs as png in folder.
}


# EDA----

outfile <- paste0("IPM_EDA.html")
outDir1 <- paste0("output/outputEDA", i)
outDir2 <- paste0("output/outputEDA_ez", i)
#if(!dir.exists(paste0("output/", outDir)))dir.create(outDir)
rmarkdown::render("IPM_EDA.Rmd", output_file = outfile, output_dir = outDir1) # makes a dashboard
