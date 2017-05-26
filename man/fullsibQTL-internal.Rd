\name{fullsibQTL-internal}

\alias{calc_probs}

\alias{calc_genoprobs}

\alias{create.map}

\alias{colin_check}

\alias{geno_class_identify}

\alias{model_select_lg}

\alias{model_select}

\alias{prob2cof}

\alias{unlink2cof}

\alias{get_phase}

\alias{draw.map2}

\title{Internal fullsibQTL functions}

\description{
  Functions from \pkg{fullsibQTL} not to be directly called by the user.
}

\usage{
calc_probs(fullsib, step, error.prob, map.function, condIndex)

calc_genoprobs(dist, geno, type, phase, step, error.prob, map.function)

create.map(dist,phase,type,step)

colin_check(probs_1pos, condIndex)

geno_class_identify(colin.number)

model_select_lg(fullsib, ws, LG)

model_select(fullsib, ws)

prob2cof(prob, colin.type)

unlink2cof(mk.unmap, obj.name, pheno.index)

get_phase(model, alpha.p, lod.ap, alpha.q, lod.aq, probs)

draw.map2(map.list, horizontal = FALSE, names = FALSE, grid = FALSE,
cex.mrk = 1, cex.grp = 0.75, trait) 

}

\author{
  Rodrigo Gazaffi, all functions except:
  
  Marcelo Mollinari and Rodrigo Gazaffi \code{calc_probs} and
  \code{calc_genoprobs}, \code{draw.map2}

  Gabriel R. A. Margarido and Rodrigo Gazaffi \code{colin_check}

}

\seealso{
  \code{fullsibQTL},

  \code{onemap}
}

\keyword{internal}

