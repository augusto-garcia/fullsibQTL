######################################################################
#
# zzz.R
#
# copyright (c) 2001, Rodrigo Gazaffi
# written March, 2011
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the onemap-QTL package
#
# .First.lib is run when the package is loaded with library(onemap-QTL)
#
######################################################################

.First.lib <- function(lib, pkg)
  library.dynam("outbredQTL", pkg, lib)

# end of zzz.R

