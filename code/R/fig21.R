rm(list = ls())

# INTRO ------------------------------------------------------------------------

# __ list of required packages ====
toInstall <- c("extrafont")


# __ verify and install uninstalled packages ====
for (pkg in toInstall) {
    if (!require(pkg, character.only = TRUE)) {
        install.packages(pkg, repos = "http://cran.us.r-project.org")
    }
}

# __ load Verdana font ====
if (!("Verdana" %in% fonts())) {
    font_import(prompt=FALSE)
    loadfonts()
}

plot.offline = FALSE
if (file.exists("OFFLINE")) { plot.offline = TRUE }

#Clean console
cat('\f')       #Clean console
graphics.off()  #Close all previously opened windows

source("derivatives.r")

allcomponents()
