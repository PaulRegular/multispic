
## Running into issues when compiling on build.
## Process hangs on ld.exe; have to close it using Task Manager.
## TMB::compile works:
wd <- getwd()
setwd(file.path(wd, "src"))
TMB::compile("multispic.cpp", framework = "TMBad")
setwd(wd)
rm(wd)

## Run above and Install without preclean

# TMB::compile ("multispic.cpp", "-O1 -g", DLLFLAGS = "")
# TMB :: gdbsource ("debug_tmb.R" , interactive = TRUE )
