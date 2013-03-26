package edu.gsu.cs.kgem.model


/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/17/13
 * Time: 8:24 PM
 * To change this template use File | Settings | File Templates.
 */
class EM(gens: List[Genotype], reads: List[Read]) {
  val rs = (0 until reads.size)
  val gs = (0 until gens.size)
  val eps = 0.005

  var h_rs = initHrs
  var freqs = Array.fill(gens.size) {
    1.0 / gens.size
  }
  var rFreqs = {
    val d = new Array[Double](reads.size)
    val s = (reads.map(r => r.freq) sum)
    reads.zipWithIndex foreach (e => {
      d(e._2) = e._1.freq / s
    })
    d
  }

  def initHrs: Array[Array[Double]] = {
    val hrs = Array.fill[Double](gens.size, reads.size) {
      1.0
    }
    if (gens.size == 0 || reads.size == 0) {
      return hrs
    }
    var g = 0

    while (g < gens.size) {
      val gen = gens(g)
      val lm = gen.data.length
      var i = 0
      while (i < lm) {
        for (c <- KGEM.table(i)) {
          val multiplier = gen.data(i)(c._1)
          for (r <- c._2) {
            val j = r._2
            hrs(g)(j) *= multiplier
          }
        }
        i += 1
      }
      g += 1
    }
    hrs
  }

  def run = {
    while (mStep(eStep) > eps) {}
    gs foreach (e => gens(e).freq = freqs(e))
  }

  def eStep = {
    val pqrs = Array.ofDim[Double](gens.size, reads.size)
    val rSums = new Array[Double](reads.size)
    rs foreach (r => {
      gs foreach (g => {
        pqrs(g)(r) = freqs(g) * h_rs(g)(r)
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

  def mStep(pqrs: Array[Array[Double]]): Double = {
    val nfqs = Array.tabulate[Double](gens.size)((i) => {
      rs.map(k => pqrs(i)(k) * rFreqs(k)).sum
    })
    val change = ((for (i <- 0 until gens.size) yield Math.abs(nfqs(i) - freqs(i))) max)
    for (i <- gs) {
      freqs(i) = nfqs(i)
    }
    println(change)
    change
  }
}
