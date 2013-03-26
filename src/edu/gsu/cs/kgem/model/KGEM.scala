package edu.gsu.cs.kgem.model

import collection.mutable

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/23/13
 * Time: 8:26 PM
 * To change this template use File | Settings | File Templates.
 */
object KGEM {
  private var reads = List[Read]()
  private var zreads = reads.zipWithIndex
  private var em = new EM(List[Genotype](), List[Read]())
  var table = new mutable.MutableList[Map[String, List[(Read, Int)]]]()

  def initReads(reads: List[Read]) = {
    this.reads = reads
    zreads = reads.zipWithIndex
    val ml = reads.map(r => r.end).max
    var i = 0
    val nucls = Genotype.sMap.keys
    while (i < ml) {
      table += ((for (nucl <- nucls)
      yield (nucl, zreads.filter(r => r._1.seq(i).equals(nucl(0))))).toMap)
      i += 1
    }
  }

  def run(gens: List[Genotype]) = {
    var collupse = gens
    var collupsed = gens.size
    do {
      collupsed = collupse.size
      runKgem(collupse)
      val col = collupse.map(g => (g.toIntegralString, g)).toMap
      collupsed -= col.size
      collupse = col.values.toList
      println("KGEM collupsed " + collupsed + " genotypes")
    } while (collupsed > 0)
    collupse
  }

  def runKgem(gens: List[Genotype]) = {
    for (g <- gens) g.convergen = false

    while (!gens.forall(g => g.convergen)) {
      val st = System.currentTimeMillis
      rounding(gens)
      runEM(gens)
      alleleFreqEstimation(gens)
      println("KGEM iteration done in " + ((System.currentTimeMillis - st) * 1.0 / 60000) + " minutes")
    }
  }

  private def rounding(gens: List[Genotype]) = {
    for (g <- gens) g.round
  }

  private def runEM(gens: List[Genotype]) = {
    em = new EM(gens, reads)
    em.run
  }


  private def alleleFreqEstimation(gens: List[Genotype]) = {
    val pqrs = em.eStep
    for (g <- gens.zipWithIndex) doAlleleFreqEstimation(g._1, pqrs(g._2))
  }

  private def doAlleleFreqEstimation(g: Genotype, pqs: Array[Double]): Unit = {
    if (g.convergen) return
    val prev = g.toIntegralString
    for (e <- g.data.zipWithIndex) {
      for (v <- e._1.keys) e._1(v) = 0
      val idx = e._2
      val rs = table(idx)
      for (n <- rs.keys) {
        val res = rs(n)
        for (r <- res)
          e._1(n) += pqs(r._2)
      }
    }
    g.convergen = prev.equals(g.toIntegralString)
  }
}
