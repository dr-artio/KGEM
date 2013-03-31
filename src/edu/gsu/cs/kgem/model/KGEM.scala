package edu.gsu.cs.kgem.model

import collection.mutable
import util.Random

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
  private var tr = 0.0005
  var table = new mutable.MutableList[Map[String, List[(Read, Int)]]]()

  def initSeeds(n: Int): List[Genotype] = {
    val seeds = sample(reads, n)
    return seeds.map(s => new Genotype(s.seq)).toList
  }

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
    var collapse = gens
    var collapsed = gens.size
    do {
      collapsed = collapse.size
      runKgem(collapse)
      val col = collapse.map(g => (g.toIntegralString, g)).toMap
      collapse = col.values.toList
      collapse = thresholdClean(collapse, tr)
      collapsed -= collapse.size
      println("KGEM collapsed %d genotypes".format(collapsed))
    } while (collapsed > 0)
    collapse
  }

  def initThreshold(tr: Double) {
    this.tr = tr
  }

  private def thresholdClean(gens: List[Genotype], tr: Double) = {
    val cleaned = gens.filter(g => g.freq >= tr)
    cleaned
  }

  private def runKgem(gens: List[Genotype]) = {
    for (g <- gens) g.convergen = false

    while (!gens.forall(g => g.convergen)) {
      val st = System.currentTimeMillis
      rounding(gens)
      runEM(gens)
      alleleFreqEstimation(gens)
      println("KGEM iteration done in %.2f minutes".format(((System.currentTimeMillis - st) * 1.0 / 60000)))
    }
  }

  private def rounding(gens: List[Genotype]) = {
    for (g <- gens) g.round
  }

  def runEM(gens: List[Genotype]) = {
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

  /**
   * Method for choosing random sample of size @size from the
   * collection of objects. Returns the whole list is @size
   * is greater than size of the original collection.
   * @param iter
   * Any iterable collection
   * @param size
   * Size of the required sample
   * @tparam T
   * Generic parameter (Class of objects in list)
   * @return
   * List of randomly chosen elements from collection
   * of specified size @size
   */
  private def sample[T](iter: Iterable[T], size: Int) = {
    if (iter.size < size) iter.toList
    var res = new mutable.MutableList[T]()
    val rnd = new Random(System.currentTimeMillis)
    var needed = size
    var len = iter.size
    val iterator = iter.iterator
    while (needed > 0 && iterator.hasNext) {
      val item = iterator.next
      if (rnd.nextInt(len) < needed) {
        res += item
        needed -= 1
      }
      len -= 1
    }
    res
  }

}
