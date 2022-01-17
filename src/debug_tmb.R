
wd <- getwd()
setwd(file.path(wd, "src"))
TMB::compile("multispic.cpp")
setwd(wd)
rm(wd)

# TMB::compile ("multispic.cpp", "-O1 -g", DLLFLAGS = "")
# TMB :: gdbsource ("debug_tmb.R" , interactive = TRUE )
