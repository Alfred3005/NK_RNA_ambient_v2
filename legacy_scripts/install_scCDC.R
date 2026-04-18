# 🧬 PHOENIX_PROTOCOL: scCDC R-Installation Fix
# Safely install scCDC into a personal library.

# 1. Define local library path
personal_lib <- "~/R/library"
if (!dir.exists(personal_lib)) {
    dir.create(personal_lib, recursive = TRUE)
}

# 2. Add to library paths
.libPaths(c(personal_lib, .libPaths()))

# 3. Check for devtools
if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "http://cran.us.r-project.org", lib = personal_lib)
}

# 4. Install scCDC
message("--- Installing scCDC from GitHub into ", personal_lib, " ---")
options(unzip = "internal") # Avoid some system-level unzip issues
devtools::install_github(
    "ZJU-UoE-CCW-LAB/scCDC", 
    upgrade = "never", 
    lib = personal_lib,
    force = TRUE
)

message("--- scCDC Installation Successful! ---")
