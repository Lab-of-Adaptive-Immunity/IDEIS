########################################################################
#
# License:
#
# IDEIS (c) by Lab of Adaptive Immunity from Institute of Molecular Genetics of the Czech Academy of Sciences
#
# IDEIS and all its parts are licensed under a
# Creative Commons Attribution 4.0 International License.
#
# You should have received a copy of the license along with this
# work. If not, see <https://creativecommons.org/licenses/by/4.0/>.
#
########################################################################
#
# Acknowledgements:
#
# This project was supported by the National Institute of Virology and 
# Bacteriology (Programme EXCELES, LX22NPO5103 to Ondrej Stepanek) - 
# funded by the European Union - Next Generation EU.
#
########################################################################
#
# Description:
#
# This is a small script to extract whitelist from .rds file to use it with IDEIS.
#
########################################################################

library(Seurat)
library(magrittr)

args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2 || length(args) > 2) {
  stop('Incorrect number of arguments detected (whitelist reading). A path to input directory, to output directory and sequencing specs are needed.')
}

input.file.path <- args[1]
output.file.path <- args[2]

rds.file <- readRDS(input.file.path)
write(gsub('-1', '', colnames(rds.file)), output.file.path)
