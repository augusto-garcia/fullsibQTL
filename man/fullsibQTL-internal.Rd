\name{fullsibQTL-internal}

\alias{calc.probs}

\alias{calc.genoprobs}

\alias{create.map}

\alias{colin.check}

\alias{geno.class.identify}

\alias{model.select.lg}

\alias{model.select}

\alias{prob2cof}

\alias{unlink2cof}

\alias{get.phase}

\alias{draw.map2}

\title{Internal fullsibQTL functions}

\description{
  Functions from \pkg{fullsibQTL} not to be directly called by the user.
}

\usage{
calc.probs(fullsib, step, error.prob, map.function, condIndex)

calc.genoprobs(dist, geno, type, phase, step, error.prob, map.function)

create.map(dist,phase,type,step)

colin.check(probs_1pos, condIndex)

geno.class.identify(colin.number)

model.select.lg(fullsib, ws, LG)

model.select(fullsib, ws)

prob2cof(prob, colin.type)

unlink2cof(mk.unmap, obj.name, pheno.index)

get.phase(model, alpha.p, lod.ap, alpha.q, lod.aq, probs)

draw.map2(map.list, horizontal = FALSE, names = FALSE, grid = FALSE,
cex.mrk = 1, cex.grp = 0.75, trait) 

}

\author{
  Rodrigo Gazaffi, all functions except:
  
  Marcelo Mollinari and Rodrigo Gazaffi \code{calc.probs} and
  \code{calc.genoprobs}, \code{draw.map2}

  Gabriel R. A. Margarido and Rodrigo Gazaffi \code{colin.check}

}

\seealso{
  \code{fullsibQTL},

  \code{onemap}
}

\keyword{internal}

