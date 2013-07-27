package edu.gsu.cs.kgem.model

import collection.mutable
import edu.gsu.cs.kgem.model.estimation.{EMMAP, EM}

/**
 * Created with IntelliJ IDEA.
 * User: aartyomenko
 * Date: 3/23/13
 * Time: 8:26 PM
 * To change this template use File | Settings | File Templates.
 */
object KGEM {
  private var reads: List[Read] = List[Read]()
  private var zreads = reads.zipWithIndex
  private var em = new EM(List[Genotype](), List[Read]())
  private var tr = 0.0005
  var table = new mutable.MutableList[Map[String, Iterable[(Read, Int)]]]()
  var loglikelihood = 0.0
  var maximumAP = 0.0

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

  def run(gens: Iterable[Genotype], alpha: Double) = {
    var collapse = gens
    var collapsed = gens.size
    var genMap = new mutable.HashMap[Genotype, Int]()
    do {
      collapsed = collapse.size
      runKgem(collapse, alpha)
      val col = collapse.map(g => (g.toIntegralString, g)).toMap
      collapse = col.values.toList
      collapse = thresholdClean(collapse, tr)
      collapse.zipWithIndex.foreach(pair => genMap += pair)
      collapsed -= collapse.size
      println("KGEM collapsed %d genotypes".format(collapsed))
    } while (collapsed > 0)
    calcLogLikelihood(collapse, genMap, alpha)
    collapse
  }

  def initThreshold(tr: Double) {
    this.tr = tr
  }

  private def calcLogLikelihood(gens: Iterable[Genotype], genIdxMap: mutable.HashMap[Genotype, Int], alpha: Double) = {
    val pqrs = em.eStep
    var ll = 0.0
    for (readpair <- zreads) {
      val read = readpair._1
      val idx = readpair._2
      ll += (read.freq * math.log(gens.map(gen => gen.freq * pqrs(genIdxMap(gen))(idx)).sum))
    }
    if (alpha == 0.0) {
      this.loglikelihood = ll
    } else {
      this.maximumAP = ll
    }
  }

  private def thresholdClean(gens: Iterable[Genotype], tr: Double) = {
    val cleaned = gens.filter(g => g.freq >= tr)
    cleaned
  }

  private def runKgem(gens: Iterable[Genotype], alpha: Double) = {
    for (g <- gens) g.convergen = false
    var i = 1
    while (!gens.forall(g => g.convergen) && i <= 5) {
      val st = System.currentTimeMillis
      rounding(gens)
      runEM(gens, alpha)
      alleleFreqEstimation(gens)
      println("KGEM iteration #%d done in %.2f minutes".format(i, ((System.currentTimeMillis - st) * 1.0 / 60000)))
      i += 1
    }
  }

  private def rounding(gens: Iterable[Genotype]) = {
    for (g <- gens) g.round
  }

  def runEM(gens: Iterable[Genotype], alpha: Double) = alpha match {
    case 0.0 => {
      em = new EM(gens.toList, reads)
      em.run
    }
    case _ => {
      em = new EMMAP(gens.toList, reads.toList, alpha)
      em.run
    }
  }


  private def alleleFreqEstimation(gens: Iterable[Genotype]) = {
    val pqrs = em.eStep
    for (g <- gens.zipWithIndex.par) doAlleleFreqEstimation(g._1, pqrs(g._2))
  }

  private def doAlleleFreqEstimation(g: Genotype, pqs: Array[Double]): Unit = {
    if (g.convergen) return
    val prev = g.toIntegralString
    for (e <- g.data.zipWithIndex.par) {
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
