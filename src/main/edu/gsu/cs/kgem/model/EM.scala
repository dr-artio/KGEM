package edu.gsu.cs.kgem.model

import edu.gsu.cs.kgem.exec.log

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/17/13
 * Time: 8:24 PM
 * Object for performing Expectation Maximization (EM)
 * estimation for prescribed edu.gsu.cs.kgem.model.
 */
object EM {
  var eps = 0.001
}

class EM(gens: List[Genotype], reads: List[Read]) {
  val rs = (0 until reads.size).view
  val gs = (0 until gens.size).view
  val eps = EM.eps

  var h_rs = initHrs
  var freqs = Array.fill(gens.size) {
    1.0 / gens.size
  }
  val total = reads.map(r => r.freq).sum
  var rFreqs = {
    val d = new Array[Double](reads.size)
    val s = total
    reads.zipWithIndex foreach (e => {
      d(e._2) = e._1.freq / s
    })
    d
  }

  /**
   * Initialize h_rs in two dimensional grid of
   * values.
   * @return
   * Two dim array with h_rs
   */
  private def initHrs: Array[Array[Double]] = {
    val hrs = Array.fill[Double](gens.size, reads.size) {
      1.0
    }
    if (gens.size == 0 || reads.size == 0) {
      return hrs
    }
    var g = 0

    for (gen <- gens.par) {
      val g = gens.indexOf(gen)
      val lm = gen.data.length
      var i = 0
      while (i < lm) {
        for (c <- KGEM.table(i)) {
          if (gen.data(i).contains(c._1.charAt(0))) {
            val multiplier = gen.data(i)(c._1.charAt(0))
            for (r <- c._2) {
              val j = r._2
              hrs(g)(j) *= multiplier
            }
          }
        }
        i += 1
      }
    }
    log("h_rs initialized.")
    hrs
  }

  /**
   * Main method for performing EM algorithm
   * Call after initialization.
   */
  def run = {
    while (mStep(eStep) > eps) {}
    gs foreach (e => gens(e).freq = freqs(e))
  }

  /**
   * Expectation step.
   * For each pair q and r p_qr=f_q*h_qr / (sum_qi_r(f_qi*h_qir))
   * @return
   * p_qrs probabilities of emitting read by corresponding
   * genotype
   */
  def eStep = {
    val pqrs = Array.ofDim[Double](gens.size, reads.size)
    val rSums = new Array[Double](reads.size)
    rs foreach (r => {
      gs foreach (g => {
        var h = h_rs(g)(r)
        if (h == Double.PositiveInfinity) h = Double.MaxValue
        pqrs(g)(r) = freqs(g) * h
        rSums(r) += pqrs(g)(r)
      })
    })
    rs foreach (r => {
      val s = rSums(r)
      if (s > 0) gs foreach (g => {
        pqrs(g)(r) /= s
      })
    })
    pqrs
  }

  /**
   * Maximization step.
   * For each qsps f_q=sum_q_ri(p_q_ri*f_ri) / sum_rj(f_rj)
   * @param pqrs
   * p_qrs probabilities of emitting read by corresponding
   * genotype
   * @return
   * Max change of frequency.
   */
  def mStep(pqrs: Array[Array[Double]]): Double = {
    val nfqs = Array.tabulate[Double](gens.size)((i) => {
      rs.view.map(k => pqrs(i)(k) * rFreqs(k)).sum
    })
    val sum = nfqs.sum
    if (sum != 1.0) gs.foreach(x => nfqs(x) /= sum)
    val change = gs.map(x => Math.abs(nfqs(x) - freqs(x))).max
    for (i <- gs) {
      freqs(i) = nfqs(i)
    }
    log("%f".format(change))
    change
  }
}
