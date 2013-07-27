package edu.gsu.cs.kgem.model.estimation

import cc.mallet.util.Dirichlet
import edu.gsu.cs.kgem.model.{Read, Genotype}

/**
 * Created with IntelliJ IDEA.
 * User: nicholas
 * Date: 7/17/13
 * Time: 12:24 PM
 * To change this template use File | Settings | File Templates.
 */
/**
 *  User: nicholas
 *  Date: 7/17/13
 *  Time: 12:24 PM
 */
class EMMAP(gens: List[Genotype], reads: List[Read], alpha: Double) extends EM(gens, reads) {
  private val diri = new Dirichlet(gens.size, alpha)
  freqs = diri.nextDistribution()

  override def mStep(pqrs: Array[Array[Double]]): Double = {
    val nfqs = Array.tabulate[Double](gens.size)((i) => {
      (rs.map(k => pqrs(i)(k) * reads(k).freq).sum + this.alpha) / (total + (reads.size * this.alpha))
    })
    val change = ((for (i <- 0 until gens.size) yield Math.abs(nfqs(i) - freqs(i))) max)
    for (i <- gs) {
      freqs(i) = nfqs(i)
    }
    println("%f".format(change))
    change
  }
}
